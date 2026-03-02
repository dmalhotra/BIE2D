#ifndef _SCTL_ODE_SOLVER_TXX_
#define _SCTL_ODE_SOLVER_TXX_

#include <stdio.h>                   // for printf
#include <algorithm>                 // for max, min

#include "sctl/common.hpp"           // for Long, Integer, SCTL_ASSERT, SCTL...
#include "sctl/ode-solver.hpp"       // for SDC
#include "sctl/comm.hpp"             // for Comm (ptr only), CommOp
#include "sctl/comm.txx"             // for Comm::Allreduce, Comm::Comm
#include "sctl/iterator.hpp"         // for Iterator, ConstIterator
#include "sctl/iterator.txx"         // for Iterator::operator[]
#include "sctl/lagrange-interp.hpp"  // for LagrangeInterp
#include "sctl/lagrange-interp.txx"  // for LagrangeInterp::Interpolate
#include "sctl/math_utils.hpp"       // for QuadReal, cos, operator*, operator-
#include "sctl/math_utils.txx"       // for const_pi, cos, pow, machine_eps
#include "sctl/matrix.hpp"           // for Matrix
#include "sctl/matrix.txx"           // for Matrix::operator[], Matrix::Matr...
#include "sctl/quadrule.hpp"         // for ChebQuadRule
#include "sctl/quadrule.txx"         // for ChebQuadRule::ComputeNdsWts
#include "sctl/static-array.hpp"     // for StaticArray
#include "sctl/vector.hpp"           // for Vector
#include "sctl/vector.txx"           // for Vector::Vector<ValueType>, Vecto...

namespace sctl {

  template <class Real> void SDC<Real>::test_one_step(const Integer Order) {
    auto ref_sol = [](Real t) { return cos<Real>(-t); };
    auto fn = [](Vector<Real>* dudt, const Vector<Real>& u) {
      (*dudt)[0] = -u[1];
      (*dudt)[1] = u[0];
    };

    const SDC<Real> ode_solver(Order);
    Real t = 0.0, dt = 1.0e-1;
    Vector<Real> u, u0(2);
    u0[0] = 1.0;
    u0[1] = 0.0;
    while (t < 10.0) {
      Real error_interp, error_picard;
      ode_solver(&u, dt, u0, fn, 0.0, &error_interp, &error_picard);
      { // Accept solution
        u0 = u;
        t = t + dt;
      }

      printf("t = %e;  ", t);
      printf("u = %e;  ", u0[0]);
      printf("error = %e;  ", ref_sol(t) - u0[0]);
      printf("time_step_error_estimate = %e;  \n", std::max(error_interp, error_picard));
    }
  }

  template <class Real> void SDC<Real>::test_adaptive_solve(const Integer Order, const Real tol) {
    auto ref_sol = [](Real t) { return cos(-t); };
    auto fn = [](Vector<Real>* dudt, const Vector<Real>& u) {
      (*dudt)[0] = -u[1];
      (*dudt)[1] = u[0];
    };

    Vector<Real> u, u0(2);
    u0[0] = 1.0; u0[1] = 0.0;
    Real T = 10.0, dt = 1.0e-1;

    SDC<Real> ode_solver(Order);
    Real t = ode_solver.AdaptiveSolve(&u, dt, T, u0, fn, tol, nullptr, true);

    if (t == T) {
      printf("u = %e;  ", u[0]);
      printf("error = %e;  \n", ref_sol(T) - u[0]);
    }
  }

  template <class Real> SDC<Real>::SDC(const Integer Order_, const Comm& comm_, const bool parallel_picard_) : order(Order_), comm(comm_), parallel_picard(parallel_picard_) {
    SCTL_ASSERT(order >= 2); // TODO: use explicit Euler if order == 1
    SetPicardIter(2*order, order);

    #ifdef SCTL_QUAD_T
    using ValueType = QuadReal;
    #else
    using ValueType = long double;
    #endif

    auto second_kind_cheb_nds = [](const Integer Order) {
      Vector<ValueType> x_cheb(Order);
      for (Long i = 0; i < Order; i++) {
        x_cheb[i] = 0.5 - 0.5 * cos(const_pi<ValueType>() * i / (Order - 1));
      }
      return x_cheb;
    };
    const auto nds0 = second_kind_cheb_nds(order); // TODO: use Gauss-Lobatto nodes
    SCTL_ASSERT(nds0.Dim() == order);

    const auto build_trunc_matrix = [](const Vector<ValueType>& nds0, const Vector<ValueType>& nds1) {
      const Integer order = nds0.Dim();
      const Integer TRUNC_Order = nds1.Dim();

      Matrix<ValueType> Minterp0(order, TRUNC_Order);
      Matrix<ValueType> Minterp1(TRUNC_Order, order);
      Vector<ValueType> interp0(order*TRUNC_Order, Minterp0.begin(), false);
      Vector<ValueType> interp1(TRUNC_Order*order, Minterp1.begin(), false);
      LagrangeInterp<ValueType>::Interpolate(interp0, nds0, nds1);
      LagrangeInterp<ValueType>::Interpolate(interp1, nds1, nds0);
      Matrix<ValueType> M_error_ = (Minterp0 * Minterp1).Transpose();
      for (Long i = 0; i < order; i++) M_error_[i][i] -= 1;

      Matrix<Real> M_error(order, order);
      for (Long i = 0; i < order*order; i++) M_error[0][i] = (Real)M_error_[0][i];
      return M_error;
    };
    { // Set M_error
      Integer TRUNC_Order = order;
      if (order >= 2) TRUNC_Order = order - 1;
      if (order >= 6) TRUNC_Order = order - 1;
      if (order >= 9) TRUNC_Order = order - 1;
      M_error = build_trunc_matrix(nds0, second_kind_cheb_nds(TRUNC_Order));
      M_error_half = build_trunc_matrix(nds0, second_kind_cheb_nds(order/2));
    }
    { // Set M_time_step
      Vector<ValueType> qx, qw;
      ChebQuadRule<ValueType>::ComputeNdsWts(&qx, &qw, order);
      const Matrix<ValueType> Mw(order, 1, (Iterator<ValueType>)qw.begin(), false);
      SCTL_ASSERT(qw.Dim() == order);
      SCTL_ASSERT(qx.Dim() == order);

      Matrix<ValueType> Minterp(order, order), M_time_step_(order, order);
      Vector<ValueType> interp(order*order, Minterp.begin(), false);
      for (Integer i = 0; i < order; i++) {
        LagrangeInterp<ValueType>::Interpolate(interp, nds0, qx*nds0[i]);
        Matrix<ValueType> M_time_step_i(order,1, M_time_step_[i], false);
        M_time_step_i = Minterp * Mw * nds0[i];
      }

      M_time_step.ReInit(order, order);
      for (Long i = 0; i < order*order; i++) M_time_step[0][i] = (Real)M_time_step_[0][i];
    }
    { // Set nds
      nds.ReInit(order);
      for (Long i = 0; i < order; i++) {
        nds[i] = (Real)nds0[i];
      }
    }
  }

  template <class Real> Integer SDC<Real>::Order() const { return order; }

  // solve u = u0 + \int_0^{dt} F(u)
  template <class Real> void SDC<Real>::operator()(Vector<Real>* u, const Real dt, const Matrix<Real>& Mu0, const FnBatch& F, const Real tol_picard, Real* error_interp, Real* error_picard, Integer* iter_count, Matrix<Real>* delta_u_substep) const {
    bool verbose = true;

    const Long DOF = Mu0.Dim(1);
    SCTL_ASSERT(Mu0.Dim(0) == 1 || Mu0.Dim(0) == order);

    const Integer Nbuff = 1000;
    StaticArray<Real,Nbuff> buff;
    StaticArray<Integer,50> failed_flag_buff;
    SCTL_ASSERT(order<50);

    Matrix<Real> Mu;
    Matrix<Real> Mf0, Mf1;
    Matrix<Real> Mv, Mv_change;
    Vector<Real> picard_err;
    if (Nbuff < 1*order*DOF) Mu.ReInit(order, DOF);
    else Mu.ReInit(order, DOF, buff + 0*order*DOF, false);
    if (Nbuff < 2*order*DOF) Mf0.ReInit(order, DOF);
    else Mf0.ReInit(order, DOF, buff + 1*order*DOF, false);
    if (Nbuff < 3*order*DOF) Mf1.ReInit(order, DOF);
    else Mf1.ReInit(order, DOF, buff + 2*order*DOF, false);
    if (Nbuff < 4*order*DOF) Mv.ReInit(order, DOF);
    else Mv.ReInit(order, DOF, buff + 3*order*DOF, false);
    if (Nbuff < 5*order*DOF) Mv_change.ReInit(order, DOF);
    else Mv_change.ReInit(order, DOF, buff + 4*order*DOF, false);
    if (Nbuff < 5*order*DOF+max_picard_iter) picard_err.ReInit(max_picard_iter);
    else picard_err.ReInit(max_picard_iter, buff + 5*order*DOF, false);

    { // Evaluate Mf0 at Mu0
      const Long batch_size = Mu0.Dim(0);
      Vector<Integer> failed_flag(batch_size, failed_flag_buff, false);
      Matrix<Real> Mf0_(batch_size, DOF, Mf0.begin(), false);
      F(&Mf0_, &failed_flag, Mu0, 0, 0);
      SCTL_ASSERT(!failed_flag[0]); // should not fail at initial condition

      for (Long i = 1; i < order; i++) {
        if (i >= failed_flag.Dim() || failed_flag[i]) { // use previous solution
          for (Long j = 0; j < DOF; j++) {
            Mf0[i][j] = Mf0[i-1][j];
          }
        }
      }
    }
    for (Long j = 0; j < DOF; j++) { // Set Mu(0,:), Mf1(0,:)
      Mu[0][j] = Mu0[0][j];
      Mf1[0][j] = Mf0[0][j];
    }
    Matrix<Real>::GEMM(Mv, M_time_step, Mf0);

    Vector<Real> r(DOF);
    Long failed_flag_ = 0, picard_iter = 0;
    for (picard_iter = 0; picard_iter < max_picard_iter; picard_iter++) { // Picard iteration
      failed_flag_ = 0;
      if (parallel_picard) {
        for (Long i = 1; i < order; i++) { // Mu <-- u0 + Mv * dt
          for (Long j = 0; j < DOF; j++) {
            Mu[i][j] = Mu0[0][j] + Mv[i][j] * dt;
          }
        }

        const Long batch_size = order-1;
        Matrix<Real> Mf1_(batch_size, DOF, Mf1[1], false);
        Vector<Integer> failed_flag(batch_size, failed_flag_buff, false);
        F(&Mf1_, &failed_flag, Matrix<Real>(batch_size,DOF,Mu[1],false), picard_iter, 1);
        for (Long i = 1; i < order; i++) { // fall back to Mf0 if evaluation fails a some sub-step
          if (failed_flag[i-1]) {
            failed_flag_ = i;
            for (Long j = 0; j < DOF; j++) {
              Mf1[i][j] = Mf0[i][j];
            }
          }
        }
      } else {
        r = 0;
        StaticArray<Integer,1> failed_flag_buff;
        Vector<Integer> failed_flag(1, failed_flag_buff, false);
        for (Long i = 1; i < order; i++) { // correction sub-steps
          const Vector<Real> f0_0(DOF, Mf0[i-1], false);
          const Vector<Real> f1_0(DOF, Mf1[i-1], false);
          Vector<Real> v_1(DOF, Mv[i], false);
          Matrix<Real> u_1(1, DOF, Mu[i], false);
          Matrix<Real> f1_1(1, DOF, Mf1[i], false);

          for (Long j = 0; j < DOF; j++) {
            r[j] += (f1_0[j] - f0_0[j]) * (nds[i]-nds[i-1]); // forward-Euler correction
            u_1[0][j] = Mu0[0][j] + (v_1[j] + r[j]) * dt;
          }

          F(&f1_1, &failed_flag, u_1, picard_iter, i);
          if (failed_flag[0]) { // use previous solution if F evaluation failed
            failed_flag_ = i;
            for (Long j = 0; j < DOF; j++) {
              Mf1[i][j] = Mf0[i][j];
            }
          }
        }
      }
      Mf0.Swap(Mf1);

      Mv_change.Swap(Mv);
      Matrix<Real>::GEMM(Mv, M_time_step, Mf0);
      Mv_change -= Mv;
      picard_err[picard_iter] = max_norm(Mv_change) * dt;
      if (verbose && !comm.Rank()) std::cout<<"SDC: picard_iter = " << picard_iter << ";  picard_err = " << picard_err[picard_iter] << "; delta_u = " << max_norm(Mv) * dt << '\n';
      if (!failed_flag_ && picard_err[picard_iter] <= tol_picard) break; // converged
      if (picard_iter-picard_stagnate_steps>0 && picard_err[picard_iter] > picard_err[picard_iter-picard_stagnate_steps]) break; // stagnated
    }

    if (verbose && !comm.Rank() && (failed_flag_ || picard_err[picard_iter] > tol_picard)) { // did not converge
      if (failed_flag_) std::cout<<"SDC: Picard iteration evaluation failed at sub-step " << failed_flag_ << ".\n";
      if (picard_iter == max_picard_iter) std::cout<<"SDC: Picard iteration did not converge in " << max_picard_iter << " iterations.\n";
      else std::cout<<"SDC: Picard iteration stagnated at picard_iter = " << picard_iter << ".\n";
      std::cout<<"SDC: Picard iteration error history: ";
      for (Long i = 0; i < std::min(picard_iter+1,max_picard_iter); i++) std::cout<<picard_err[i]<<"  ";
      std::cout<<'\n';
    }

    if (u) { // Set output u <-- u0 + Mv[order-1] * dt
      if (failed_flag_) {
        u->ReInit(0);
      } else {
        if (u->Dim() != DOF) u->ReInit(DOF);
        for (Long j = 0; j < DOF; j++) {
          (*u)[j] = Mu0[0][j] + Mv[order-1][j] * dt;
        }
      }
    }

    if (error_picard != nullptr) {
      (*error_picard) = picard_err[std::min<Long>(picard_iter, max_picard_iter-1)];
    }
    if (error_interp != nullptr) {
      Matrix<Real>& err = Mv_change; // reuse memory
      Matrix<Real>::GEMM(err, M_error, Mv); // truncation error of interpolant coefficients
      (*error_interp) = max_norm(err) * dt;
    }
    if (iter_count != nullptr) {
      (*iter_count) = 1+std::min<Long>(picard_iter, max_picard_iter-1);
    }
    if (delta_u_substep != nullptr) {
      if (delta_u_substep->Dim(0) != order || delta_u_substep->Dim(1) != DOF) delta_u_substep->ReInit(order, DOF);
      for (Long i = 0; i < order; i++) {
        for (Long j = 0; j < DOF; j++) {
          (*delta_u_substep)[i][j] = Mv[i][j] * dt;
        }
      }
    }
  }

  template <class Real> void SDC<Real>::operator()(Vector<Real>* u, const Real dt, const Vector<Real>& u0, const Fn0& F, const Real tol_picard, Real* error_interp, Real* error_picard, Integer* iter_count, Matrix<Real>* u_substep) const {
    const auto fn = [&F](Matrix<Real>* dudt, Vector<Integer>* failed_flag, const Matrix<Real>& u, const Integer correction_idx, const Integer substep_offset) {
      const Long N_substep = u.Dim(0), N = u.Dim(1);
      if (dudt->Dim(0) != N_substep || dudt->Dim(1) != N) dudt->ReInit(N_substep, N);
      if (failed_flag && failed_flag->Dim() != N_substep) failed_flag->ReInit(N_substep);

      for (Long i = 0; i < N_substep; i++) {
        const Vector<Real> u_(N, (Iterator<Real>)u[i], false);
        Vector<Real> dudt_(N, (*dudt)[i], false);
        F(&dudt_, u_, correction_idx, substep_offset+i);
        if (failed_flag) (*failed_flag)[i] = (dudt_.Dim() == 0); // mark failed evaluations with empty output
      }
    };
    this->operator()(u, dt, Matrix<Real>(1,u0.Dim(),(Iterator<Real>)u0.begin(),false), fn, tol_picard, error_interp, error_picard, iter_count, u_substep);
  }

  template <class Real> void SDC<Real>::operator()(Vector<Real>* u, const Real dt, const Vector<Real>& u0, const Fn1& F, const Real tol_picard, Real* error_interp, Real* error_picard, Integer* iter_count, Matrix<Real>* u_substep) const {
    const auto fn = [&F](Vector<Real>* dudt, const Vector<Real>& u, const Integer correction_idx, const Integer substep_idx) {
      F(dudt, u);
    };
    this->operator()(u, dt, u0, fn, tol_picard, error_interp, error_picard, iter_count, u_substep);
  }

  template <class Real> Real SDC<Real>::AdaptiveSolve(Vector<Real>* u, Real dt, const Real T, const Vector<Real>& u0, const FnBatch& F, const Real tol, const MonitorFn* monitor_callback, bool continue_with_errors, Real* error) const {
    const Real eps = machine_eps<Real>();
    const auto nds0 = GetNodes();
    const Long DOF = u0.Dim();

    Real t = 0;
    Real error_ = 0;
    Vector<Real> u_;
    Matrix<Real> delta_u_substep;
    Matrix<Real> u0_(1, DOF, (Iterator<Real>)u0.begin());
    while (t < T && dt > eps*T) {
      bool accept = false;
      Integer picard_iter;
      Real error_interp, error_picard;
      Real tol_ = std::max<Real>(tol/T, (tol-error_)/(T-t));
      (*this)(&u_, dt, u0_, F, tol_*dt, &error_interp, &error_picard, &picard_iter, &delta_u_substep);

      const Real max_norm_u = max_norm(u_);
      const Real error_interp_half = (u_.Dim() ? pow<Real>(max_norm(M_error_half*delta_u_substep)/max_norm_u, (Real)1.8) * max_norm_u : 0);

      if (!comm.Rank()) std::cout<<"Adaptive time-step: " << std::scientific << std::setw(10) <<t<<' '<<dt<<' '<<picard_iter<<" "<<error_interp/dt<<" "<<error_interp_half/dt<<' '<<error_picard/dt<<' '<<error_/tol<<'\n';
      if (u_.Dim() && (error_interp < tol_*dt || error_interp_half < error_interp) && (error_picard < tol_*dt)) { // Accept solution
        // u_.Dim()                           // SDC time-step succeeded
        // && (
        //   error_interp < tol_*dt           // interpolant error tolerance reached
        //   ||
        //   error_interp_half < error_interp // interpolant coefficients stagnated (this usually cannot be fixed by reducing time-step size, so we continue)
        // ) && (
        //   error_picard < tol_*dt           // picard-iteration error tolerance reached
        // )

        t = t + dt;
        accept = true;
        error_ += std::max<Real>(error_interp, error_picard);
        SCTL_ASSERT_MSG(continue_with_errors || error_ < tol, "Could not solve ODE to the requested tolerance.");
        if (monitor_callback) (*monitor_callback)(t, dt, u_);
      }

      const Real dt_picard = [this,&u_,&dt,&error_picard,&tol_,&picard_iter](){
        if (!u_.Dim()) return dt*(Real)0.5; // aborted
        if (picard_iter < max_picard_iter && error_picard > tol_*dt) return dt*(Real)0.5; // diverged
        return dt * (max_picard_iter*0.66)/picard_iter * log(error_picard)/log(tol_*dt); // adjust step size aiming for picard_iter to be around 66% of max_picard_iter
      }();
      const Real dt_interp = (error_interp_half < error_interp ?
                              std::min<Real>(1.5*dt, 0.9*dt * pow<Real>(error_interp/error_interp_half, 1/(Real)(order))) : // adjust time-step size to match stagnation error
                              std::max<Real>(0.5*dt, 0.9*dt * pow<Real>((tol_*dt)   /error_interp_half, 1/(Real)(order)))); // Adjust time-step size (Quaife, Biros - JCP 2016)
      if (!comm.Rank()) std::cout<<"current dt="<<dt<<"  new dt_picard="<<dt_picard<<"  new dt_interp="<<dt_interp<<'\n'; ///////////////////////////////
      const Real dt_new = std::min<Real>(T-t, std::min(dt_interp, dt_picard));

      { // Build u0_ for next step
        if (accept) { // extrapolate from sub-step solutions
          //const Vector<Real> nds1{1};
          //const Vector<Real> nds1{0, 1};
          //const Vector<Real> nds1{0, 0.5, 1};
          //const Vector<Real> nds1{0, 1/(Real)3, 2/(Real)3, 1};
          //const Vector<Real> nds1{0, 0.25, 0.5, 0.75, 1};
          const Vector<Real> nds1{0, 0.2, 0.4, 0.6, 0.8, 1};

          if (nds1.Dim() == 1) {
            u0_ = Matrix<Real>(1,DOF,u_.begin());
          } else {
            static const auto Minterp0_t = [this,&nds0,&nds1](){
              Matrix<Real> Minterp0(order, nds1.Dim());
              Vector<Real> interp0(order*nds1.Dim(), Minterp0.begin(), false);
              LagrangeInterp<Real>::Interpolate(interp0, nds0, nds1);
              return Minterp0.Transpose();
            }();

            const Vector<Real> nds2 = 1 + nds0 * (dt_new/dt);
            Matrix<Real> Minterp1(nds1.Dim(), order);
            Vector<Real> interp1(nds1.Dim()*order, Minterp1.begin(), false);
            LagrangeInterp<Real>::Interpolate(interp1, nds1, nds2);

            Vector<Real> u0(DOF, u0_.begin());
            u0_ = Minterp1.Transpose() * (Minterp0_t * delta_u_substep);
            for (Long i = 0; i < order; i++) {
              for (Long j = 0; j < DOF; j++) {
                u0_[i][j] += u0[j];
              }
            }
          }
        } else { // interpolate from sub-step solutions
          const Vector<Real> nds1 = nds0 * (dt_new/dt);
          Matrix<Real> Minterp(order, order);
          Vector<Real> interp(order*order, Minterp.begin(), false);
          LagrangeInterp<Real>::Interpolate(interp, nds0, nds1);

          Vector<Real> u0(DOF, u0_.begin());
          u0_ = Minterp.Transpose() * delta_u_substep;
          for (Long i = 0; i < order; i++) {
            for (Long j = 0; j < DOF; j++) {
              u0_[i][j] += u0[j];
            }
          }
        }
      }
      dt = dt_new;
    }
    if (t < T || error_ > tol) SCTL_WARN("Could not solve ODE to the requested tolerance.");
    if (error != nullptr) (*error) = error_;

    (*u) = Vector<Real>(DOF, u0_.begin(), false);
    return t;
  }

  template <class Real> Real SDC<Real>::AdaptiveSolve(Vector<Real>* u, Real dt, const Real T, const Vector<Real>& u0, const Fn0& F, Real tol, const MonitorFn* monitor_callback, bool continue_with_errors, Real* error) const {
    const auto fn = [&F](Matrix<Real>* dudt, Vector<Integer>* failed_flag, const Matrix<Real>& u, const Integer correction_idx, const Integer substep_offset) {
      const Long N_substep = u.Dim(0), N = u.Dim(1);
      if (dudt->Dim(0) != N_substep || dudt->Dim(1) != N) dudt->ReInit(N_substep, N);
      if (failed_flag && failed_flag->Dim() != N_substep) failed_flag->ReInit(N_substep);

      for (Long i = 0; i < N_substep; i++) {
        const Vector<Real> u_(N, (Iterator<Real>)u[i], false);
        Vector<Real> dudt_(N, (*dudt)[i], false);
        F(&dudt_, u_, correction_idx, substep_offset+i);
        if (failed_flag) (*failed_flag)[i] = (dudt_.Dim() == 0); // mark failed evaluations with empty output
      }
    };
    return AdaptiveSolve(u, dt, T, u0, fn, tol, monitor_callback, continue_with_errors, error);
  }

  template <class Real> Real SDC<Real>::AdaptiveSolve(Vector<Real>* u, Real dt, const Real T, const Vector<Real>& u0, const Fn1& F, Real tol, const MonitorFn* monitor_callback, bool continue_with_errors, Real* error) const {
    const auto fn = [&F](Vector<Real>* dudt, const Vector<Real>& u, const Integer correction_idx, const Integer substep_idx) {
      F(dudt, u);
    };
    return AdaptiveSolve(u, dt, T, u0, fn, tol, monitor_callback, continue_with_errors, error);
  }

  template <class Real> template <class Container> Real SDC<Real>::max_norm(const Container& M) const {
    StaticArray<Real,2> max_val{0,0};
    for (const auto x : M) max_val[0] = std::max<Real>(max_val[0], fabs((Real)x));
    comm.Allreduce((ConstIterator<Real>)max_val, (Iterator<Real>)max_val+1, 1, CommOp::MAX);
    return max_val[1];
  }
}

#endif // _SCTL_ODE_SOLVER_TXX_
