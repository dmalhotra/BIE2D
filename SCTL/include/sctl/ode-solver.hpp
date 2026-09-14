#ifndef _SCTL_ODE_SOLVER_HPP_
#define _SCTL_ODE_SOLVER_HPP_

#include <functional>       // for function

#include "sctl/common.hpp"  // for Integer, sctl
#include "sctl/comm.hpp"    // for Comm
#include "sctl/comm.txx"    // for Comm::Self
#include "sctl/matrix.hpp"  // for Matrix
#include "sctl/vector.hpp"  // for Vector

namespace sctl {

/**
 * Implements spectral deferred correction (SDC) solver for ordinary differential equations (ODEs).
 */
template <class Real> class SDC {
  public:

    /// The function type to specify the RHS of the ODE.
    using FnBatch = std::function<void(Matrix<Real>* dudt, Vector<Integer>* failed_flag, const Matrix<Real>& u, const Integer correction_idx, const Integer substep_offset)>;

    /// The function type to specify the RHS of the ODE.
    using Fn0 = std::function<void(Vector<Real>* dudt, const Vector<Real>& u, const Integer correction_idx, const Integer substep_idx)>;

    /// The function type to specify the RHS of the ODE.
    using Fn1 = std::function<void(Vector<Real>* dudt, const Vector<Real>& u)>;

    /// Callback function type.
    using MonitorFn = std::function<void(Real t, Real dt, const Vector<Real>& u)>;

    /**
     * Constructor
     *
     * @param[in] order the order of the method.
     *
     * @param[in] comm the communicator.
     *
     * @param[in] parallel_picard whether to parallelize the Picard iterations across sub-steps.
     */
    explicit SDC(const Integer order, const Comm& comm = Comm::Self(), const bool parallel_picard = true);

    void SetPicardIter(Integer max_iter, Integer stagnate_steps) {
      max_picard_iter = max_iter;
      picard_stagnate_steps = stagnate_steps;
    }

    const Vector<Real>& GetNodes() const { return nds; }

    /**
     * @return order of the method.
     */
    Integer Order() const;

    /**
     * Apply one step of spectral deferred correction (SDC).
     * Compute: \f$ u = u_0 + \int_0^{dt} F(u) \f$
     *
     * @param[out] u the solution
     * @param[in] dt the step size
     *
     * @param[in] u0 matrix of size (M, DOF), where M is either 1 or order.
     * u0(0,:) is the initial value and u0(i,:) for i>0 is an optional initial
     * guess for the solution at the i-th sub-step.
     *
     * @param[in] F the function du/dt
     * @param[in] tol_picard the tolerance for stopping Picard iterations
     * @param[out] error_interp an estimate of the truncation error of the solution interpolant
     * @param[out] error_picard the Picard iteration error
     * @param[out] iter_count number of Picard iterations
     * @param[out] delta_u_substep the value of (u(t)-u0) at each sub-step
     */
    void operator()(Vector<Real>* u, const Real dt, const Matrix<Real>& u0, const FnBatch& F, const Real tol_picard = 0, Real* error_interp = nullptr, Real* error_picard = nullptr, Integer* iter_count = nullptr, Matrix<Real>* delta_u_substep = nullptr) const;

    void operator()(Vector<Real>* u, const Real dt, const Vector<Real>& u0, const Fn0& F, const Real tol_picard = 0, Real* error_interp = nullptr, Real* error_picard = nullptr, Integer* iter_count = nullptr, Matrix<Real>* delta_u_substep = nullptr) const;

    /**
     * Apply one step of spectral deferred correction (SDC).
     * Compute: \f$ u = u_0 + \int_0^{dt} F(u) \f$
     *
     * @param[out] u the solution
     * @param[in] dt the step size
     *
     * @param[in] u0 matrix of size (M, DOF), where M is either 1 or order.
     * u0(0,:) is the initial value and u0(i,:) for i>0 is an optional initial
     * guess for the solution at the i-th sub-step.
     *
     * @param[in] F the function du/dt
     * @param[in] tol_picard the tolerance for stopping Picard iterations
     * @param[out] error_interp an estimate of the truncation error of the solution interpolant
     * @param[out] error_picard the Picard iteration error
     * @param[out] iter_count number of Picard iterations on exit (or -1 if terminated)
     * @param[out] u_substep the solution at each substep
     */
    void operator()(Vector<Real>* u, const Real dt, const Vector<Real>& u0, const Fn1& F, const Real tol_picard = 0, Real* error_interp = nullptr, Real* error_picard = nullptr, Integer* iter_count = nullptr, Matrix<Real>* u_substep = nullptr) const;

    /**
     * Solve ODE adaptively to required tolerance.
     * Compute: \f$ u = u_0 + \int_0^{T} F(u) \f$
     *
     * @param[out] u the final solution
     * @param[in] dt the initial step size guess
     * @param[in] T the final time
     * @param[in] u0 the initial value
     * @param[in] F the function du/dt
     * @param[in] tol the required solution tolerance
     * @param[in] monitor_callback a callback function called after each accepted time-step
     * @param[in] continue_with_errors tries to compute the best solution even if the required tolerance cannot be satisfied.
     * @param[out] error estimate of the final output error
     *
     * @return the final time (should equal T if no errors)
     */
    Real AdaptiveSolve(Vector<Real>* u, Real dt, const Real T, const Vector<Real>& u0, const FnBatch& F, Real tol, const MonitorFn* monitor_callback = nullptr, bool continue_with_errors = false, Real* error = nullptr, bool adaptive_step_size = true) const;

    Real AdaptiveSolve(Vector<Real>* u, Real dt, const Real T, const Vector<Real>& u0, const Fn0& F, Real tol, const MonitorFn* monitor_callback = nullptr, bool continue_with_errors = false, Real* error = nullptr) const;

    /**
     * Solve ODE adaptively to required tolerance.
     * Compute: \f$ u = u_0 + \int_0^{T} F(u) \f$
     *
     * @param[out] u the final solution
     * @param[in] dt the initial step size guess
     * @param[in] T the final time
     * @param[in] u0 the initial value
     * @param[in] F the function du/dt
     * @param[in] tol the required solution tolerance
     * @param[in] monitor_callback a callback function called after each accepted time-step
     * @param[in] continue_with_errors tries to compute the best solution even if the required tolerance cannot be satisfied.
     * @param[out] error estimate of the final output error
     *
     * @return the final time (should equal T if no errors)
     */
    Real AdaptiveSolve(Vector<Real>* u, Real dt, const Real T, const Vector<Real>& u0, const Fn1& F, Real tol, const MonitorFn* monitor_callback = nullptr, bool continue_with_errors = false, Real* error = nullptr) const;

    /**
     * This is an example for how to use the SDC class.
     */
    static void test_one_step(const Integer Order = 5);

    /**
     * This example shows adaptive time-stepping with the SDC class.
     */
    static void test_adaptive_solve(const Integer Order = 5, const Real tol = 1e-5);

  private:

    template <class Container> Real max_norm(const Container& M) const;

    Matrix<Real> M_time_step, M_error, M_error_half;
    Vector<Real> nds;
    Integer order;
    Comm comm;

    Integer max_picard_iter;
    Integer picard_stagnate_steps;
    bool parallel_picard;
};

}

#endif // _SCTL_ODE_SOLVER_HPP_
