#include <bie2d.hpp>
using namespace sctl;

const char data_path[]="data2/";
constexpr Integer ElemOrder = 64;
constexpr Integer MatInterpOrder = 32;

template <class ValueType, class Real> Matrix<ValueType> convert_mat(Matrix<Real> M) {
  Matrix<ValueType> M_;
  M_.ReInit(M.Dim(0), M.Dim(1));
  for (long i = 0; i < M.Dim(0)*M.Dim(1); i++) M_[0][i] = (ValueType)M[0][i];
  return M_;
};

template <class Real> Real cond(const Matrix<Real>& M) {
  //if (M.Dim(0)*M.Dim(1) > 1e6) return -1;
  if (!M.Dim(0) || !M.Dim(1)) return -1;
  Matrix<Real> U, S, Vt;
  Matrix<Real>(M).SVD(U, S, Vt);
  Real Smax = S[0][0], Smin = S[0][0];
  const Long N = std::min(S.Dim(0), S.Dim(1));
  for (Long i = 0; i < N; i++) {
    Smax = std::max<Real>(Smax, fabs(S[i][i]));
    Smin = std::min<Real>(Smin, fabs(S[i][i]));
  }
  return Smax / Smin;
}

template <class Real> Matrix<Real> pinv(const Matrix<Real>& M, const Real tol) {
  using ValueType = Real;
  Matrix<ValueType> M_(M.Dim(0), M.Dim(1));
  for (Long i = 0; i < M.Dim(0)*M.Dim(1); i++) M_[0][i] = (ValueType)M[0][i];
  auto Minv_ = M_.pinv((ValueType)tol);

  Matrix<Real> Minv(Minv_.Dim(0), Minv_.Dim(1));
  for (Long i = 0; i < Minv.Dim(0)*Minv.Dim(1); i++) Minv[0][i] = (Real)Minv_[0][i];
  return Minv;
}

template <class Real> void print_sing(const Matrix<Real>& M) {
  Matrix<Real> U, S, Vt;
  Matrix<Real>(M).SVD(U,S,Vt);
  for (Long i = 0; i < std::min(S.Dim(0),S.Dim(1)); i++) {
    std::cout<<S[i][i]<<'\n';
  }
  std::cout<<'\n';
}

template <class Real> Matrix<Real> apply_inv(const Matrix<Real>& A, const Matrix<Real>& B, const Real tol) {
  Matrix<Real> U, S, Vt;
  Matrix<Real>(B).SVD(U, S, Vt);

  Real S_max = 0;
  SCTL_ASSERT(S.Dim(0) == S.Dim(1));
  for (Long i = 0; i < S.Dim(0); i++) {
    S_max = std::max<Real>(S_max, fabs(S[i][i]));
  }
  for (Long i = 0; i < S.Dim(0); i++) {
    if (fabs(S[i][i]) > S_max*tol) {
      S[i][i] = 1/S[i][i];
    } else {
      S[i][i] = 0;
    }
  }

  return ((A * Vt.Transpose()) * S) * U.Transpose();
}


template <class Real, class Disc> Vector<Real> SolveCapacitance(const Vector<Disc>& disc_lst, const bool use_gmres = true, const bool use_sqrtscal = true, const Real tol = machine_eps<Real>()*64, const Long max_iter = 1000) {
  using KerSL = Laplace2D_FxU;
  using KerDL = Laplace2D_DxU;

  Vector<Real> Xt;
  GetGeom(disc_lst, &Xt);

  Matrix<Real> Msl, Mdl;
  LayerPotentialMatrix<KerSL>(Msl, disc_lst, Xt, tol);
  LayerPotentialMatrix<KerDL>(Mdl, disc_lst, Xt, tol);

  Vector<Real> sqrt_w, invsqrt_w;
  { // Set sqrt_w, insqrt_w
    BoundaryIntegralWts(sqrt_w, disc_lst);
    if (use_sqrtscal) {
      for (auto& w : sqrt_w) w = sqrt<Real>(w);
    } else {
      sqrt_w = 1;
    }
    invsqrt_w = 1 / sqrt_w;
  }

  Matrix<Real> M_bie;
  { // Set M_bie
    M_bie = Mdl + Msl;
    for (Long i = 0; i < M_bie.Dim(0); i++) M_bie[i][i] += 0.5;
    for (Long i = 0; i < M_bie.Dim(0); i++) {
      for (Long j = 0; j < M_bie.Dim(1); j++) {
        M_bie[i][j] *= invsqrt_w[i] * sqrt_w[j];
      }
    }
  }
  auto BIOp = [&M_bie](Vector<Real>* Ax_, const Vector<Real>& x_) {
    const Long N = x_.Dim();
    SCTL_ASSERT(Ax_->Dim() == N);

    const Matrix<Real> x(1, N, (Iterator<Real>)x_.begin(), false);
    Matrix<Real> Ax(1, N, Ax_->begin(), false);
    Matrix<Real>::GEMM(Ax, x, M_bie);
  };

  // Setup boundary integral operator: K = 1/2 + D + S (exterior limit)
  auto BIOp_ = [&Msl,&Mdl,&sqrt_w,&invsqrt_w](Vector<Real>* Ax, const Vector<Real>& x_) {
    const auto x = x_ * invsqrt_w;

    Vector<Real> Usl(x.Dim()), Udl(x.Dim()); Usl = 0; Udl = 0;
    Matrix<Real> Usl_(1, Usl.Dim(), Usl.begin(), false);
    Matrix<Real> Udl_(1, Udl.Dim(), Udl.begin(), false);
    Matrix<Real>::GEMM(Usl_, Matrix<Real>(1, x.Dim(), (Iterator<Real>)x.begin(), false), Msl);
    Matrix<Real>::GEMM(Udl_, Matrix<Real>(1, x.Dim(), (Iterator<Real>)x.begin(), false), Mdl);
    //LayerPotential<KerSL>(Usl, disc_lst, Xt, x, tol);
    //LayerPotential<KerDL>(Udl, disc_lst, Xt, x, tol);
    (*Ax) = 0.5*x + Udl + Usl;
    (*Ax) *= sqrt_w;
  };

  Vector<Real> U0; // Dirichlet boundary conditions
  for (Long i = 0; i < disc_lst.Dim(); i++) { // Set U0
    const Long N = disc_lst[i].NodeCount();
    for (Long j = 0; j < N; j++) {
      U0.PushBack(i+1); // U0 = i for i-th disc
    }
  }

  Vector<Real> sigma; // unknown density
  if (use_gmres) { // GMRES
    SCTL_UNUSED(BIOp_);
    ParallelSolver<Real> solver; // Linear solver (GMRES)
    solver(&sigma, BIOp, U0*sqrt_w, tol, max_iter);
    sigma *= invsqrt_w;
  } else { // Direct
    std::cout<< "Condition number  = " << cond(M_bie) << '\n';
    Matrix<Real> Minv_bie = pinv(M_bie, tol);
    for (Long i = 0; i < M_bie.Dim(0); i++) {
      for (Long j = 0; j < M_bie.Dim(1); j++) {
        Minv_bie[i][j] *= sqrt_w[i] * invsqrt_w[j];
      }
    }

    const Matrix<Real> U0_(1, U0.Dim(), U0.begin(), false);
    Matrix<Real> sigma_ = U0_ * Minv_bie;
    sigma.ReInit(sigma_.Dim(1), sigma_.begin());
  }

  // Return charge on each disk
  Vector<Real> Q;
  Long offset = 0;
  for (auto& disc : disc_lst) {
    const Long N = disc.NodeCount();
    Vector<Real> sigma_(N, sigma.begin() + offset, false);

    Vector<Real> q;
    disc.BoundaryIntegralDirect(q, sigma_);
    Q.PushBack(q[0]);

    offset += N;
  }
  return Q;
}

template <class Real, class Disc> void SolveElastance(const Vector<Disc>& disc_lst, const Vector<Real> Q, const bool use_gmres = true, const bool use_sqrtscal = true, const Real tol = machine_eps<Real>()*64, const Long max_iter = 1000) {
  using KerSL = Laplace2D_FxU;
  using KerDL = Laplace2D_DxU;

  Vector<Real> Xt;
  GetGeom(disc_lst, &Xt);

  Matrix<Real> Mdl;
  LayerPotentialMatrix<KerDL>(Mdl, disc_lst, Xt, tol);

  Vector<Real> sqrt_w, invsqrt_w;
  { // Set sqrt_w, insqrt_w
    BoundaryIntegralWts(sqrt_w, disc_lst);
    if (use_sqrtscal) {
      for (auto& w : sqrt_w) w = sqrt<Real>(w);
    } else {
      sqrt_w = 1;
    }
    invsqrt_w = 1 / sqrt_w;
  }

  Matrix<Real> M_bie;
  { // Set M_bie
    M_bie = Mdl;
    for (Long i = 0; i < M_bie.Dim(0); i++) M_bie[i][i] += 0.5;

    Long offset = 0;
    for (auto& disc : disc_lst) {
      Vector<Real> W;
      disc.BoundaryIntegralWts(W);
      for (Long i = 0; i < W.Dim(); i++) {
        for (Long j = 0; j < W.Dim(); j++) {
          M_bie[offset+i][offset+j] += W[i];
        }
      }
      offset += W.Dim();
    }

    for (Long i = 0; i < M_bie.Dim(0); i++) {
      for (Long j = 0; j < M_bie.Dim(1); j++) {
        M_bie[i][j] *= invsqrt_w[i] * sqrt_w[j];
      }
    }
  }
  auto BIOp = [&M_bie](Vector<Real>* Ax_, const Vector<Real>& x_) {
    const Long N = x_.Dim();
    SCTL_ASSERT(Ax_->Dim() == N);

    const Matrix<Real> x(1, N, (Iterator<Real>)x_.begin(), false);
    Matrix<Real> Ax(1, N, Ax_->begin(), false);
    Matrix<Real>::GEMM(Ax, x, M_bie);
  };

  // Setup boundary integral operator: K = 1/2 + D (exterior limit)
  auto BIOp_ = [&disc_lst,&Mdl,&sqrt_w,&invsqrt_w](Vector<Real>* Ax, const Vector<Real>& x_) {
    const auto x = x_ * invsqrt_w;

    Vector<Real> Udl(x.Dim()); Udl = 0;
    //LayerPotential<KerDL>(Udl, disc_lst, Xt, x, tol);
    Matrix<Real> Udl_(1, Udl.Dim(), Udl.begin(), false);
    Matrix<Real>::GEMM(Udl_, Matrix<Real>(1, x.Dim(), (Iterator<Real>)x.begin(), false), Mdl);
    (*Ax) = 0.5*x + Udl;

    Long offset = 0;
    for (auto& disc : disc_lst) {
      const Long N = disc.NodeCount();
      const Vector<Real> x_(N, (Iterator<Real>)x.begin() + offset, false);
      Vector<Real> Ax_(N, Ax->begin() + offset, false);

      Vector<Real> x0;
      disc.BoundaryIntegralDirect(x0, x_);
      Ax_ += x0[0];

      offset += N;
    }
    (*Ax) *= sqrt_w;
  };

  Vector<Real> nu, S_nu;
  SCTL_ASSERT(disc_lst.Dim() == Q.Dim());
  for (Long i = 0; i < disc_lst.Dim(); i++) {
    Real nu0 = Q[i] / (2*const_pi<Real>()*disc_lst[i].Radius());
    for (Long j = 0; j < disc_lst[i].NodeCount(); j++) {
      nu.PushBack(nu0);
    }
  }
  LayerPotential<KerSL>(S_nu, disc_lst, Xt, nu, tol); // Snu <-- S[nu]

  Vector<Real> sigma; // unknown density
  if (use_gmres) { // GMRES
    SCTL_UNUSED(BIOp_);
    ParallelSolver<Real> solver; // Linear solver (GMRES)
    solver(&sigma, BIOp, S_nu*(-1)*sqrt_w, tol, max_iter);
    sigma *= invsqrt_w;
  } else { // Direct
    std::cout<< "Condition number  = " << cond(M_bie) << '\n';
    Matrix<Real> Minv_bie = pinv(M_bie, tol);
    for (Long i = 0; i < M_bie.Dim(0); i++) {
      for (Long j = 0; j < M_bie.Dim(1); j++) {
        Minv_bie[i][j] *= sqrt_w[i] * invsqrt_w[j];
      }
    }

    const Matrix<Real> S_nu_(1, S_nu.Dim(), S_nu.begin(), false);
    Matrix<Real> sigma_ = (S_nu_*(-1)) * Minv_bie;
    sigma.ReInit(sigma_.Dim(1), sigma_.begin());
  }

  Vector<Real> D_sigma;
  LayerPotential<KerDL>(D_sigma, disc_lst, Xt, sigma, tol); // D_sigma <-- D[sigma]
  const Vector<Real> U = 0.5*sigma + D_sigma + S_nu;

  // print error
  Long offset = 0;
  Real max_err = 0;
  for (Long i = 0; i < disc_lst.Dim(); i++) {
    const Long N = disc_lst[i].NodeCount();
    for (Long j = 0; j < N; j++) {
      max_err = std::max<Real>(max_err, fabs(U[offset+j]-i-1));
    }
    offset += N;
  }
  std::cout<<"Error = "<<max_err<<'\n';
}








template <class Real, class Disc> void CapaElasOpMatrix(Matrix<Real>& Kcapa, Matrix<Real>& Kelas, const Vector<Disc>& disc_lst, const Real tol, const bool use_sqrtscal) {
  using KerSL = Laplace2D_FxU;
  using KerDL = Laplace2D_DxU;

  Vector<Real> Xt;
  GetGeom(disc_lst, &Xt);

  Matrix<Real> Msl, Mdl;
  LayerPotentialMatrix<KerSL>(Msl, disc_lst, Xt, tol);
  LayerPotentialMatrix<KerDL>(Mdl, disc_lst, Xt, tol);

  Vector<Real> sqrt_w, invsqrt_w;
  { // Set sqrt_w, insqrt_w
    BoundaryIntegralWts(sqrt_w, disc_lst);
    if (use_sqrtscal) {
      for (auto& w : sqrt_w) w = sqrt<Real>(w);
    } else {
      sqrt_w = 1;
    }
    invsqrt_w = 1 / sqrt_w;
  }

  { // Set Kcapa
    Kcapa = Mdl + Msl;
    for (Long i = 0; i < Kcapa.Dim(0); i++) Kcapa[i][i] += 0.5;
    for (Long i = 0; i < Kcapa.Dim(0); i++) {
      for (Long j = 0; j < Kcapa.Dim(1); j++) {
        Kcapa[i][j] *= invsqrt_w[i] * sqrt_w[j];
      }
    }
  }
  { // Set Kelas
    Kelas = Mdl;
    for (Long i = 0; i < Kelas.Dim(0); i++) Kelas[i][i] += 0.5;

    Long offset = 0;
    for (const auto& disc : disc_lst) {
      Vector<Real> W;
      disc.BoundaryIntegralWts(W);
      for (Long i = 0; i < W.Dim(); i++) {
        for (Long j = 0; j < W.Dim(); j++) {
          Kelas[offset+i][offset+j] += W[i];
        }
      }
      offset += W.Dim();
    }

    for (Long i = 0; i < Kelas.Dim(0); i++) {
      for (Long j = 0; j < Kelas.Dim(1); j++) {
        Kelas[i][j] *= invsqrt_w[i] * sqrt_w[j];
      }
    }
  }
}

template <class Real> class CompressedOperators {
  public:

    CompressedOperators(const Real R, const Real eps0, const Real eps1, const Long Order) : R_(R), eps0_(eps0), eps1_(eps1), Order_(Order) {
      const Real log_eps0 = log<Real>(eps0_);
      const Real log_eps1 = log<Real>(eps1_);
      log_nds = log_eps0 + (log_eps1-log_eps0) * ChebQuadRule<Real>::ComputeNds(Order_);

      Rcapa_lst.ReInit(Order_);
      Relas_lst.ReInit(Order_);
      for (Long i = 0; i < Order_; i++) {
        GetCompressPrecond_(Rcapa_lst[i], Relas_lst[i], exp<Real>(log_nds[i]));
      }
    }

    /**
     * Get Rcapa and Relas by interpolating using log-Chebyshev basis.
     */
    template <class ValueType> void GetCompressPrecond(Matrix<ValueType>& Rcapa, Matrix<ValueType>& Relas, const ValueType eps) {
      Vector<Real> interp_wts, trg_node(1); trg_node = log<Real>(eps);
      LagrangeInterp<Real>::Interpolate(interp_wts, log_nds, trg_node);

      Matrix<Real> Rcapa_ = Rcapa_lst[0]*0;
      Matrix<Real> Relas_ = Relas_lst[0]*0;
      for (Long i = 0; i < interp_wts.Dim(); i++) {
        Rcapa_ += Rcapa_lst[i] * interp_wts[i];
        Relas_ += Relas_lst[i] * interp_wts[i];
      }

      if (Rcapa.Dim(0)!=Rcapa_.Dim(0) || Rcapa.Dim(1)!=Rcapa_.Dim(1)) {
        Rcapa.ReInit(Rcapa_.Dim(0), Rcapa_.Dim(1));
      }
      if (Relas.Dim(0)!=Relas_.Dim(0) || Relas.Dim(1)!=Relas_.Dim(1)) {
        Relas.ReInit(Relas_.Dim(0), Relas_.Dim(1));
      }
      for (Long i = 0; i < Rcapa.Dim(0)*Rcapa.Dim(1); i++) Rcapa[0][i] = (ValueType)Rcapa_[0][i];
      for (Long i = 0; i < Relas.Dim(0)*Relas.Dim(1); i++) Relas[0][i] = (ValueType)Relas_[0][i];
    }

    /**
     * Compute Rcapa and Relas for a given R and eps.
     */
    template <class ValueType> static void CompressPrecond(Matrix<ValueType>& Rcapa, Matrix<ValueType>& Relas, Matrix<ValueType>* Kcapa_inv_, Matrix<ValueType>* Kelas_inv_, const ValueType R, const ValueType eps, const ValueType tol = (ValueType)machine_eps<Real>()) {
      const Integer adap_depth = 25; //(Integer)(log(eps)/log((Real)0.5)*2/3);
      Vector<Disc<Real,ElemOrder>> disc_lst(2);
      { // Init disc_lst
        auto theta_adap0 = [](Long depth) {
          Vector<Real> theta;
          for (Long i = 0; i <=depth; i++) theta.PushBack(0-pow<Real>(0.5, i));
          theta.PushBack(0);
          for (Long i = 0; i <=depth; i++) theta.PushBack(0+pow<Real>(0.5, depth-i));
          return theta * const_pi<Real>()/6;
        }(adap_depth);
        auto theta_adap1 = [](Long depth) {
          Vector<Real> theta;
          for (Long i = 0; i <=depth; i++) theta.PushBack(6-pow<Real>(0.5, i));
          theta.PushBack(6);
          for (Long i = 0; i <=depth; i++) theta.PushBack(6+pow<Real>(0.5, depth-i));
          return theta * const_pi<Real>()/6;
        }(adap_depth);
        disc_lst[0] = Disc<Real,ElemOrder>(-R-eps/2, 0.0, R, Vector<Real>(theta_adap0.Dim()-1, theta_adap0.begin()+0), Vector<Real>(theta_adap0.Dim()-1, theta_adap0.begin()+1));
        disc_lst[1] = Disc<Real,ElemOrder>( R+eps/2, 0.0, R, Vector<Real>(theta_adap1.Dim()-1, theta_adap1.begin()+0), Vector<Real>(theta_adap1.Dim()-1, theta_adap1.begin()+1));
      }

      const Matrix<Real> P = [](Long depth) {
        const Vector<Real> nds0 = LegQuadRule<Real>::ComputeNds(ElemOrder);
        Vector<Real> nds1(depth*ElemOrder);
        for (Long d = 0; d < depth; d++) { // Set nds1
          Vector<Real> nds1_(ElemOrder, nds1.begin() + d*ElemOrder, false);
          nds1_ = (1-pow<Real>(0.5,d)) + nds0 * pow<Real>(0.5, std::min(d+1,depth-1));
        }

        Matrix<Real> P0(nds0.Dim(), nds1.Dim());
        Vector<Real> P0_(P0.Dim(0)*P0.Dim(1), P0.begin(), false);
        LagrangeInterp<Real>::Interpolate(P0_, nds0, nds1);

        Matrix<Real> P(8*ElemOrder, 4*(depth+1)*ElemOrder); P = 0;
        for (Long i = 0; i < ElemOrder; i++) {
          P[0*ElemOrder+i][(0*(depth+1)+0)*ElemOrder+i] = 1;
          P[3*ElemOrder+i][(2*(depth+1)-1)*ElemOrder+i] = 1;
          P[4*ElemOrder+i][(2*(depth+1)+0)*ElemOrder+i] = 1;
          P[7*ElemOrder+i][(4*(depth+1)-1)*ElemOrder+i] = 1;
        }
        for (Long i = 0; i < ElemOrder; i++) {
          for (Long j = 0; j < depth*ElemOrder; j++) {
            P[1*ElemOrder+i+0][(0*(depth+1)+1)*ElemOrder+j+0] = P0[i][j];
            P[3*ElemOrder-i-1][(2*(depth+1)-1)*ElemOrder-j-1] = P0[i][j];
            P[5*ElemOrder+i+0][(2*(depth+1)+1)*ElemOrder+j+0] = P0[i][j];
            P[7*ElemOrder-i-1][(4*(depth+1)-1)*ElemOrder-j-1] = P0[i][j];
          }
        }
        return P;
      }(adap_depth);
      const Long Nc = P.Dim(0);
      const Long Nf = P.Dim(1);

      Matrix<Real> Wf, Wc_inv;
      { // Set Wf
        Vector<Real> Wf_;
        for (const auto& disc : disc_lst) {
          Vector<Real> W;
          disc.BoundaryIntegralWts(W);
          for (const auto& w : W) Wf_.PushBack(w);
        }
        SCTL_ASSERT(Wf_.Dim() == Nf);
        Wf.ReInit(Nf, Nf);
        Wf = 0;
        for (Long i = 0; i < Wf_.Dim(); i++) {
          Wf[i][i] = Wf_[i];
        }
      }
      { // Set Wc_inv
        const auto Wc = P*Wf*P.Transpose();
        Wc_inv.ReInit(Nc, Nc);
        Wc_inv = 0;
        for (Long i = 0; i < Nc; i++) {
          Wc_inv[i][i] = 1 / Wc[i][i];
        }
      }

      Matrix<Real> Kcapa, Kelas;
      CapaElasOpMatrix(Kcapa, Kelas, disc_lst, (Real)tol, false);

      auto Compress = [&P,&Wf,&Wc_inv](Matrix<ValueType>* Rcomp, Matrix<ValueType>* Kinv, const Matrix<Real>& Kfine, const Real tol) {
        Matrix<Real> U, S, Vt;
        Matrix<Real>(Kfine).SVD(U, S, Vt);

        Real S_max = 0;
        for (Long i = 0; i < S.Dim(0); i++) {
          S_max = std::max<Real>(S_max, fabs(S[i][i]));
        }
        for (Long i = 0; i < S.Dim(0); i++) {
          if (S[i][i] < tol*S_max) S[i][i] = 0;
          else S[i][i] = 1/S[i][i];
        }

        if (Kinv != nullptr) {
          Matrix<Real> Kinv_ = Vt.Transpose() * S * U.Transpose();
          if (Kinv->Dim(0) != Kinv_.Dim(0) || Kinv->Dim(1) != Kinv_.Dim(1)) {
            Kinv->ReInit(Kinv_.Dim(0), Kinv_.Dim(1));
          }
          for (Long i = 0; i < Kinv_.Dim(0)*Kinv_.Dim(1); i++) {
            (*Kinv)[0][i] = (ValueType)Kinv_[0][i];
          }
        }
        if (Rcomp != nullptr) {
          Matrix<Real> Rcomp_ = (P * Vt.Transpose()) * S * (U.Transpose() * Wf * P.Transpose() * Wc_inv);
          if (Rcomp->Dim(0) != Rcomp_.Dim(0) || Rcomp->Dim(1) != Rcomp_.Dim(1)) {
            Rcomp->ReInit(Rcomp_.Dim(0), Rcomp_.Dim(1));
          }
          for (Long i = 0; i < Rcomp_.Dim(0)*Rcomp_.Dim(1); i++) {
            (*Rcomp)[0][i] = (ValueType)Rcomp_[0][i];
          }
        }
      };
      Compress(&Rcapa, Kcapa_inv_, Kcapa, (Real)tol);
      Compress(&Relas, Kelas_inv_, Kelas, (Real)tol);

      //const Matrix<Real> Kcapa_inv = Kcapa.pinv(tol);
      //const Matrix<Real> Kelas_inv = Kelas.pinv(tol);
      //if (Kcapa_inv_ != nullptr) (*Kcapa_inv_) = Kcapa_inv;
      //if (Kelas_inv_ != nullptr) (*Kelas_inv_) = Kelas_inv;
      //Rcapa = P * Kcapa_inv * Wf * P.Transpose() * Wc_inv;
      //Relas = P * Kelas_inv * Wf * P.Transpose() * Wc_inv;
    }

    /**
     * Get Rcapa and Relas by interpolating from diadically refined Chebyshev panels.
     */
    template <class ValueType> void GetCompressPrecond_(Matrix<ValueType>& Rcapa, Matrix<ValueType>& Relas, const ValueType eps) {
      const auto cheb_nds = ChebQuadRule<Real>::ComputeNds(MatInterpOrder);
      const Long interval = (Long)(log(eps)/log((Real)0.5));

      const Real eps0 = pow<Real>(0.5, interval+0);
      const Real eps1 = pow<Real>(0.5, interval+1);
      const auto interp_nds = eps0 + (eps1-eps0)*cheb_nds;

      Vector<Matrix<Real>> Rcapa0(MatInterpOrder), Relas0(MatInterpOrder);
      for (Long i = 0; i < MatInterpOrder; i++) {
        const auto fname = PrecompFname(R_, interval, i);
        Rcapa0[i].Read((fname+"-capa").c_str());
        Relas0[i].Read((fname+"-elas").c_str());
        if (Rcapa0[i].Dim(0)*Rcapa0[i].Dim(1) == 0 || Relas0[i].Dim(0)*Relas0[i].Dim(1)==0) {
          PrecompCompressPrecond(R_, interp_nds[i], interp_nds[i], Comm::Self());
          Rcapa0[i].Read((fname+"-capa").c_str());
          Relas0[i].Read((fname+"-elas").c_str());
        }
      }

      Vector<Real> interp_wts, trg_node(1); trg_node = eps;
      LagrangeInterp<Real>::Interpolate(interp_wts, interp_nds, trg_node);

      Matrix<Real> Rcapa_ = Rcapa0[0]*0;
      Matrix<Real> Relas_ = Relas0[0]*0;
      for (Long i = 0; i < interp_wts.Dim(); i++) {
        Rcapa_ += Rcapa0[i] * interp_wts[i];
        Relas_ += Relas0[i] * interp_wts[i];
      }

      if (Rcapa.Dim(0)!=Rcapa_.Dim(0) || Rcapa.Dim(1)!=Rcapa_.Dim(1)) {
        Rcapa.ReInit(Rcapa_.Dim(0), Rcapa_.Dim(1));
      }
      if (Relas.Dim(0)!=Relas_.Dim(0) || Relas.Dim(1)!=Relas_.Dim(1)) {
        Relas.ReInit(Relas_.Dim(0), Relas_.Dim(1));
      }
      for (Long i = 0; i < Rcapa.Dim(0)*Rcapa.Dim(1); i++) Rcapa[0][i] = (ValueType)Rcapa_[0][i];
      for (Long i = 0; i < Relas.Dim(0)*Relas.Dim(1); i++) Relas[0][i] = (ValueType)Relas_[0][i];
    }

  private:

    static std::string PrecompFname(const Real R, const Long interval, const Long node_idx) {
      std::string fname = "data/Precom";
      fname += "-f" + std::to_string(sizeof(Real)*8);
      fname += "-R" + std::to_string((double)R);
      fname += "-q" + std::to_string(ElemOrder);
      fname += "-p" + std::to_string(MatInterpOrder);
      fname += "-i" + std::to_string(interval);
      fname += "-j" + std::to_string(node_idx);
      return fname;
    }

    /**
     * Compute Rcapa and Relas at nodes of diadically refined Chebyshev panels
     * and save them to file.
     */
    static void PrecompCompressPrecond(const Real R, const Real max_dist, const Real min_dist, const Comm& comm) {
      //const Real max_dist = 8e-1, min_dist = 8e-3;
      //const Real max_dist = 4e-3, min_dist = 4e-5;
      //const Real max_dist = 2e-5, min_dist = 2e-7;
      //const Real max_dist = 1e-7, min_dist = 1e-9;
      //const Real max_dist = 5e-10, min_dist = 5e-12;
      //const Real max_dist = 2.5e-12, min_dist = 2.5e-14;
      //const Real max_dist = 1.25e-14, min_dist = 1.25e-16;

      const Long np = comm.Size();
      const Long rank = comm.Rank();

      const Long idx0 = ((Long)(log(max_dist)/log((Real)0.5)+0))*MatInterpOrder;
      const Long idx1 = ((Long)(log(min_dist)/log((Real)0.5)+1))*MatInterpOrder;

      const Long idx0_ = idx0 + (idx1-idx0)*(rank+0)/np;
      const Long idx1_ = idx0 + (idx1-idx0)*(rank+1)/np;
      const auto cheb_nds = ChebQuadRule<Real>::ComputeNds(MatInterpOrder);
      for (Long i = idx0_; i < idx1_; i++) {
        const Long interval = i/MatInterpOrder;
        const Long node_idx = i%MatInterpOrder;
        const Real eps0 = pow<Real>(0.5, interval+0);
        const Real eps1 = pow<Real>(0.5, interval+1);
        const Real eps = eps0 + (eps1-eps0) * cheb_nds[node_idx];
        const auto fname = PrecompFname(R, interval, node_idx);

        Matrix<Real> Rcapa, Relas;
        Rcapa.Read((fname+"-capa").c_str());
        Relas.Read((fname+"-elas").c_str());
        if (Rcapa.Dim(0) * Rcapa.Dim(1) == 0 || Relas.Dim(0)*Relas.Dim(1)==0) {
          Matrix<Real> Kcapa_inv, Kelas_inv;
          CompressPrecond(Rcapa, Relas, &Kcapa_inv, &Kelas_inv, R, eps, machine_eps<Real>());
          Rcapa.Write((fname+"-capa").c_str());
          Relas.Write((fname+"-elas").c_str());
          //Kcapa_inv.Write((fname+"-Kinv-capa").c_str());
          //Kelas_inv.Write((fname+"-Kinv-elas").c_str());
        }
      }
    }

    Vector<Matrix<Real>> Rcapa_lst, Relas_lst;
    Vector<Real> log_nds;
    Real R_, eps0_, eps1_;
    Long Order_;
};

template <class Real> Real TestCompressPrecond(const Real R, const Real eps, const Real tol, const Long max_iter) {
  //static CompressedOperators<QuadReal> CompOp(R, atoreal<QuadReal>("1e-1"), atoreal<QuadReal>("1e-15"), 60);

  using DiscType = Disc<Real,ElemOrder>;
  Vector<DiscType> disc_lst(2);
  { // Coarse mesh
    const Long adap_depth = 1;
    auto theta_adap0 = [](Long depth) {
      Vector<Real> theta;
      for (Long i = 0; i < 10; i++) theta.PushBack(1+i);
      for (Long i = 0; i <=depth; i++) theta.PushBack(12-pow<Real>(0.5, i));
      theta.PushBack(12);
      for (Long i = 0; i <=depth; i++) theta.PushBack(12+pow<Real>(0.5, depth-i));
      return theta * const_pi<Real>()/6;
    }(adap_depth);
    auto theta_adap1 = [](Long depth) {
      Vector<Real> theta;
      for (Long i = 0; i <=depth; i++) theta.PushBack(6-pow<Real>(0.5, i));
      theta.PushBack(6);
      for (Long i = 0; i <=depth; i++) theta.PushBack(6+pow<Real>(0.5, depth-i));
      for (Long i = 0; i < 10; i++) theta.PushBack(7+(1+i));
      return theta * const_pi<Real>()/6;
    }(adap_depth);
    disc_lst[0] = DiscType(-R-eps/2, 0.0, R, Vector<Real>(theta_adap0.Dim()-1, theta_adap0.begin()+0), Vector<Real>(theta_adap0.Dim()-1, theta_adap0.begin()+1));
    disc_lst[1] = DiscType( R+eps/2, 0.0, R, Vector<Real>(theta_adap1.Dim()-1, theta_adap1.begin()+0), Vector<Real>(theta_adap1.Dim()-1, theta_adap1.begin()+1));
  }
  const Long N = NodeCount(disc_lst);

  Matrix<Real> Rcapa(N, N), Relas(N, N);
  { // Set Rcapa, Relas
    Matrix<Real> Rcapa_, Relas_;
    if (0) { // interpolate
      //CompOp.GetCompressPrecond(Rcapa_, Relas_, eps);
    } else {
      if (eps == atoreal<Real>("1e-10") && std::is_same<Real,QuadReal>::value && ElemOrder == 64) { /////////////////////////////////////
        Rcapa_.Read("Rcapa.mat");
        Relas_.Read("Relas.mat");
        if (Rcapa_.Dim(0) == 0 || Relas_.Dim(0) ==0) {
          CompressedOperators<QuadReal>::template CompressPrecond<Real>(Rcapa_, Relas_, nullptr, nullptr, R, eps);
          Rcapa_.Write("Rcapa.mat");
          Relas_.Write("Relas.mat");
        }
      } else {
      CompressedOperators<QuadReal>::template CompressPrecond<Real>(Rcapa_, Relas_, nullptr, nullptr, R, eps);
      }
    }

    Rcapa = 0;
    Relas = 0;
    for (Long i = 0; i < N; i++) {
      Rcapa[i][i] = 1;
      Relas[i][i] = 1;
    }
    for (Long i = 0; i < 8*ElemOrder; i++) {
      for (Long j = 0; j < 8*ElemOrder; j++) {
        Rcapa[N/2-4*ElemOrder+i][N/2-4*ElemOrder+j] = Rcapa_[i][j];
        Relas[N/2-4*ElemOrder+i][N/2-4*ElemOrder+j] = Relas_[i][j];
      }
    }
  }

  static Matrix<Real> Mcapa, Melas; // regularization
  auto ComputeM_identity = [&N](Matrix<Real>& Mcapa, Matrix<Real>& Melas) {
    Mcapa.ReInit(N,N);
    Melas.ReInit(N,N);
    Mcapa = 0;
    Melas = 0;
    for (Long i = 0; i < N; i++) {
      Mcapa[i][i] = 1;
      Melas[i][i] = 1;
    }
  };
  auto ComputeM_regularize = [&R,&N](Matrix<Real>& Mcapa, Matrix<Real>& Melas) { // Set Mcapa, Melas
    using ValueType = QuadReal;
    constexpr ValueType tol = machine_eps<ValueType>();
    Vector<Disc<ValueType,ElemOrder>> disc_lst(2);
    { // Init disc_lst
      auto theta_adap0 = [](Long depth) {
        Vector<ValueType> theta;
        for (Long i = 0; i <=depth; i++) theta.PushBack(0-pow<ValueType>(0.5, i));
        theta.PushBack(0);
        for (Long i = 0; i <=depth; i++) theta.PushBack(0+pow<ValueType>(0.5, depth-i));
        return theta * const_pi<ValueType>()/6;
      }(1);
      disc_lst[0] = Disc<ValueType,ElemOrder>(-R, 0.0, R, Vector<ValueType>(theta_adap0.Dim()-1, theta_adap0.begin()+0), Vector<ValueType>(theta_adap0.Dim()-1, theta_adap0.begin()+1));

      const ValueType pi_6 = const_pi<ValueType>()/6;
      const ValueType sin_pi_6 = sin<ValueType>(pi_6);
      const ValueType cos_pi_6 = sqrt<ValueType>(1-sin_pi_6*sin_pi_6);
      const ValueType RR = R*sqrt<ValueType>(sin_pi_6 * sin_pi_6 + (1-cos_pi_6) * (1-cos_pi_6));
      disc_lst[1] = Disc<ValueType,ElemOrder>(0.0, 0.0, RR, 64);

      if (1) { // adjacent panels
        Vector<ValueType> theta0, theta1;
        theta0.PushBack(-const_pi<ValueType>()/3);
        theta1.PushBack(-const_pi<ValueType>()/6);
        theta0.PushBack(const_pi<ValueType>()/6);
        theta1.PushBack(const_pi<ValueType>()/3);
        disc_lst.PushBack(Disc<ValueType,ElemOrder>(-R, 0.0, R, theta0, theta1));
      }
    }

    static const Matrix<ValueType> Minterp = []() { // Set Minterp
      Vector<ValueType> nds0 = LegQuadRule<ValueType>::ComputeNds(1*ElemOrder)*4;
      Vector<ValueType> nds1;
      { // Set nds1
        nds1.ReInit(4*ElemOrder);
        Vector<ValueType> n0(ElemOrder, nds1.begin()+0*ElemOrder, false);
        Vector<ValueType> n1(ElemOrder, nds1.begin()+1*ElemOrder, false);
        Vector<ValueType> n2(ElemOrder, nds1.begin()+2*ElemOrder, false);
        Vector<ValueType> n3(ElemOrder, nds1.begin()+3*ElemOrder, false);
        n0 = LegQuadRule<ValueType>::ComputeNds(ElemOrder) + 0;
        n1 = LegQuadRule<ValueType>::ComputeNds(ElemOrder) + 1;
        n2 = LegQuadRule<ValueType>::ComputeNds(ElemOrder) + 2;
        n3 = LegQuadRule<ValueType>::ComputeNds(ElemOrder) + 3;
      }
      const Long N0 = nds0.Dim();
      const Long N1 = nds1.Dim();

      Matrix<ValueType> Minterp0(N0, N1);
      Vector<ValueType> Minterp0_(N0*N1, Minterp0.begin(), false);
      LagrangeInterp<ValueType>::Interpolate(Minterp0_, nds0, nds1);
      return Minterp0;
    }();

    const Long N1 = 4*ElemOrder;
    const Long N2 = NodeCount(disc_lst) - N1;
    Matrix<ValueType> Mcapa__(N1, N2), Melas__(N1, N2);
    { // Set Mcapa__, Melas__
      Matrix<ValueType> Mcapa_, Melas_;
      CapaElasOpMatrix(Mcapa_, Melas_, disc_lst, tol, false);
      for (Long i = 0; i < N1; i++) {
        for (Long j = 0; j < N2; j++) {
          Mcapa__[i][j] = Mcapa_[i][N1+j];
          Melas__[i][j] = Melas_[i][N1+j];
        }
      }
    }

    //Matrix<ValueType> Q(ElemOrder, ElemOrder), S(ElemOrder, ElemOrder), Qinv;
    //{ // Set Q, Qinv
    //  const Vector<ValueType> nds0 = LegQuadRule<ValueType>::ComputeNds(1*ElemOrder);
    //  for (Long i = 0; i < ElemOrder; i++) {
    //    for (Long j = 0; j < ElemOrder; j++) {
    //      const Real t = acos<ValueType>(nds0[j]*2-1);
    //      Q[i][j] = cos<ValueType>(i*t); //*pow<ValueType>(0.5,4*i);
    //    }
    //  }
    //  Qinv = Matrix<ValueType>(Q).pinv();
    //}

    print_sing(Mcapa__);
    print_sing(Melas__);

    Mcapa__ = apply_inv<ValueType>(Mcapa__, Minterp*Mcapa__, 1e-12)*Minterp;
    Melas__ = apply_inv<ValueType>(Melas__, Minterp*Melas__, 1e-12)*Minterp;

    convert_mat<double>(Mcapa__).Write("Mcapa.mat");
    convert_mat<double>(Melas__).Write("Melas.mat");



    //Mcapa__ = apply_inv<ValueType>(Mcapa__, Minterp*Mcapa__, 1e-15) * Minterp;
    //Melas__ = apply_inv<ValueType>(Melas__, Minterp*Melas__, 1e-15) * Minterp;
    //Mcapa__ = (apply_inv<ValueType>(Mcapa__, (Q*Minterp)*Mcapa__, 1e-25) * Q) *Minterp;
    //Melas__ = (apply_inv<ValueType>(Melas__, (Q*Minterp)*Melas__, 1e-25) * Q) *Minterp;

    Mcapa.ReInit(N,N);
    Melas.ReInit(N,N);
    Mcapa = 0;
    Melas = 0;
    for (Long i = 0; i < N; i++) {
      Mcapa[i][i] = 1;
      Melas[i][i] = 1;
    }
    for (Long i = 0; i < 4*ElemOrder; i++) {
      for (Long j = 0; j < 4*ElemOrder; j++) {
        Mcapa[N/2-4*ElemOrder+i][N/2-4*ElemOrder+j] = (Real)Mcapa__[i][j];
        Mcapa[N/2-0*ElemOrder+i][N/2-0*ElemOrder+j] = (Real)Mcapa__[i][j];

        Melas[N/2-4*ElemOrder+i][N/2-4*ElemOrder+j] = (Real)Melas__[i][j];
        Melas[N/2-0*ElemOrder+i][N/2-0*ElemOrder+j] = (Real)Melas__[i][j];
      }
    }
  };
  auto ComputeM_equivalent = [&R,&eps,&N](Matrix<Real>& Mcapa, Matrix<Real>& Melas) { // Set Mcapa, Melas
    using ValueType = QuadReal;
    constexpr ValueType tol = machine_eps<ValueType>();
    Vector<Disc<ValueType,ElemOrder>> disc_lst(3);
    { // Init disc_lst
      auto theta_adap0 = [](Long depth) {
        Vector<ValueType> theta;
        for (Long i = 0; i <=depth; i++) theta.PushBack(0-pow<ValueType>(0.5, i));
        theta.PushBack(0);
        for (Long i = 0; i <=depth; i++) theta.PushBack(0+pow<ValueType>(0.5, depth-i));
        return theta * const_pi<ValueType>()/6;
      }(1);
      auto theta_adap1 = [](Long depth) {
        Vector<ValueType> theta;
        for (Long i = 0; i <=depth; i++) theta.PushBack(6-pow<ValueType>(0.5, i));
        theta.PushBack(6);
        for (Long i = 0; i <=depth; i++) theta.PushBack(6+pow<ValueType>(0.5, depth-i));
        return theta * const_pi<ValueType>()/6;
      }(1);
      disc_lst[0] = Disc<ValueType,ElemOrder>(-R-eps/2, 0.0, R, Vector<ValueType>(theta_adap0.Dim()-1, theta_adap0.begin()+0), Vector<ValueType>(theta_adap0.Dim()-1, theta_adap0.begin()+1));
      disc_lst[1] = Disc<ValueType,ElemOrder>( R+eps/2, 0.0, R, Vector<ValueType>(theta_adap1.Dim()-1, theta_adap1.begin()+0), Vector<ValueType>(theta_adap1.Dim()-1, theta_adap1.begin()+1));

      const ValueType pi_6 = const_pi<ValueType>()/6;
      const ValueType R_sin_pi_6 = R*sin<ValueType>(pi_6);
      const ValueType R_cos_pi_6 = R*cos<ValueType>(pi_6);
      const ValueType RR = sqrt<ValueType>(R_sin_pi_6 * R_sin_pi_6 + (R+eps/2-R_cos_pi_6) * (R+eps/2-R_cos_pi_6));
      disc_lst[2] = Disc<ValueType,ElemOrder>(0.0, 0.0, RR, 64);

      {
        Vector<ValueType> theta0, theta1;
        theta0.PushBack(-const_pi<ValueType>()/3);
        theta1.PushBack(-const_pi<ValueType>()/6);
        theta0.PushBack(const_pi<ValueType>()/6);
        theta1.PushBack(const_pi<ValueType>()/3);
        disc_lst.PushBack( Disc<ValueType,ElemOrder>(-R-eps/2, 0.0, R, theta0, theta1) );
      }
      {
        Vector<ValueType> theta0, theta1;
        theta0.PushBack(const_pi<ValueType>()-const_pi<ValueType>()/3);
        theta1.PushBack(const_pi<ValueType>()-const_pi<ValueType>()/6);
        theta0.PushBack(const_pi<ValueType>()+const_pi<ValueType>()/6);
        theta1.PushBack(const_pi<ValueType>()+const_pi<ValueType>()/3);
        disc_lst.PushBack( Disc<ValueType,ElemOrder>(R+eps/2, 0.0, R, theta0, theta1) );
      }
    }

    static const Matrix<ValueType> Minterp = []() { // Set Minterp
      Vector<ValueType> nds0 = LegQuadRule<ValueType>::ComputeNds(3*ElemOrder)*4;
      Vector<ValueType> nds1;
      { // Set nds1
        nds1.ReInit(4*ElemOrder);
        Vector<ValueType> n0(ElemOrder, nds1.begin()+0*ElemOrder, false);
        Vector<ValueType> n1(ElemOrder, nds1.begin()+1*ElemOrder, false);
        Vector<ValueType> n2(ElemOrder, nds1.begin()+2*ElemOrder, false);
        Vector<ValueType> n3(ElemOrder, nds1.begin()+3*ElemOrder, false);
        n0 = LegQuadRule<ValueType>::ComputeNds(ElemOrder) + 0;
        n1 = LegQuadRule<ValueType>::ComputeNds(ElemOrder) + 1;
        n2 = LegQuadRule<ValueType>::ComputeNds(ElemOrder) + 2;
        n3 = LegQuadRule<ValueType>::ComputeNds(ElemOrder) + 3;
      }
      const Long N0 = nds0.Dim();
      const Long N1 = nds1.Dim();

      Matrix<ValueType> Minterp0(N0, N1);
      Vector<ValueType> Minterp0_(N0*N1, Minterp0.begin(), false);
      LagrangeInterp<ValueType>::Interpolate(Minterp0_, nds0, nds1);

      Matrix<ValueType> Minterp(2*N0, 2*N1);
      Minterp = 0;
      for (Long i = 0; i < N0; i++) {
        for (Long j = 0; j < N1; j++) {
          Minterp[0*N0+i][0*N1+j] = Minterp0[i][j];
          Minterp[1*N0+i][1*N1+j] = Minterp0[i][j];
        }
      }
      return Minterp;
    }();

    Matrix<ValueType> Mcapa_, Melas_;
    CapaElasOpMatrix(Mcapa_, Melas_, disc_lst, tol, false);

    const Long N1 = 8*ElemOrder;
    const Long N2 = Mcapa_.Dim(1) - N1;
    Matrix<ValueType> Mcapa__(N1, N2), Melas__(N1, N2+2);
    { // Set Melas[:][N2,N2+1]
      Melas__ = 0;
      Vector<ValueType> W0, W1;
      disc_lst[0].BoundaryIntegralWts(W0);
      disc_lst[1].BoundaryIntegralWts(W1);
      for (Long i = 0; i < 4*ElemOrder; i++) {
        Melas__[0*ElemOrder+i][N2+0] = W0[i];
        Melas__[4*ElemOrder+i][N2+1] = W1[i];
      }
    }
    for (Long i = 0; i < N1; i++) {
      for (Long j = 0; j < N2; j++) {
        Mcapa__[i][j] = Mcapa_[i][N1+j];
        Melas__[i][j] = Melas_[i][N1+j];
      }
    }

    print_sing(Mcapa__);
    print_sing(Melas__);

    //Mcapa__ = apply_inv<ValueType>(Mcapa__, Mcapa__, 1e-23);
    //Melas__ = apply_inv<ValueType>(Melas__, Melas__, 1e-23);
    Mcapa__ = apply_inv<ValueType>(Mcapa__, Minterp*Mcapa__, 1e-25) * Minterp;
    Melas__ = apply_inv<ValueType>(Melas__, Minterp*Melas__, 1e-25) * Minterp;

    print_sing(Mcapa__);
    print_sing(Melas__);

    Mcapa.ReInit(N,N);
    Melas.ReInit(N,N);
    Mcapa = 0;
    Melas = 0;
    for (Long i = 0; i < N; i++) {
      Mcapa[i][i] = 1;
      Melas[i][i] = 1;
    }
    for (Long i = 0; i < 8*ElemOrder; i++) {
      for (Long j = 0; j < 8*ElemOrder; j++) {
        Mcapa[N/2-4*ElemOrder+i][N/2-4*ElemOrder+j] = (Real)Mcapa__[i][j];
        Melas[N/2-4*ElemOrder+i][N/2-4*ElemOrder+j] = (Real)Melas__[i][j];
      }
    }
  };
  if (Mcapa.Dim(0)*Mcapa.Dim(1) == 0) {
    //ComputeM_identity(Mcapa, Melas);
    ComputeM_regularize(Mcapa, Melas);
  }
  //ComputeM_equivalent(Mcapa, Melas);
  Rcapa = Rcapa * Mcapa;
  Relas = Relas * Melas;

  convert_mat<double>(Rcapa).Write("Rcapa_.mat");
  convert_mat<double>(Relas).Write("Relas_.mat");
  std::cout<<"Write R\n";

  Matrix<Real> Kcapa, Kelas; // Compressed, block-preconditioned
  { // Set Kcapa, Kelas
    Matrix<Real> Kcapa_, Kelas_;
    CapaElasOpMatrix(Kcapa_, Kelas_, disc_lst, tol, false);
    convert_mat<double>(Kcapa_).Write("Kcapa_.mat");
    convert_mat<double>(Kelas_).Write("Kelas_.mat");

    if (1) {
      //Kcapa = Rcapa * Mcapa * Kcapa_;
      //Kelas = Relas * Melas * Kelas_;
      Kcapa = Rcapa * Kcapa_;
      Kelas = Relas * Kelas_;

      for (Long i = 0; i < 8*ElemOrder; i++) {
        for (Long j = 0; j < 8*ElemOrder; j++) {
          Kcapa[N/2-4*ElemOrder+i][N/2-4*ElemOrder+j] = (i==j ? 1 : 0);
          Kelas[N/2-4*ElemOrder+i][N/2-4*ElemOrder+j] = (i==j ? 1 : 0);
        }
      }
    }

    Kcapa = Matrix<Real>(Kcapa).pinv(1e-22) * Rcapa;
    Kelas = Matrix<Real>(Kelas).pinv(1e-22) * Relas;

    auto compute_inv = [](Matrix<Real> M, const Real tol) {
      Matrix<Real> U,S,Vt;
      M.SVD(U,S,Vt);

      Real S_max = 0;
      SCTL_ASSERT(S.Dim(0) == S.Dim(1));
      for (Long i = 0; i < S.Dim(0); i++) {
        S_max = std::max<Real>(S_max, fabs(S[i][i]));
      }
      for (Long i = 0; i < S.Dim(0); i++) {
        if (fabs(S[i][i]) > S_max*tol) {
          S[i][i] = 1/S[i][i];
        } else {
          S[i][i] = 1;
        }
      }
      return (Vt.Transpose() * S) * U.Transpose();
    };
    Kcapa = compute_inv(Kcapa, 1e-22);
    Kelas = compute_inv(Kelas, 1e-22);
    //Kcapa = Matrix<Real>(Kcapa).pinv(1e-12);
    //Kelas = Matrix<Real>(Kelas).pinv(1e-12);

    { // Construct Rinv + K
    //  const auto Rcapa_inv = Matrix<Real>(Rcapa).pinv();
    //  const auto Relas_inv = Matrix<Real>(Relas).pinv();

    //  Kcapa = Rcapa_inv * Kcapa;
    //  Kelas = Relas_inv * Kelas;
    //  //Kcapa = Kcapa_;
    //  //Kelas = Kelas_;
    //  //for (Long i = 0; i < 8*ElemOrder; i++) {
    //  //  for (Long j = 0; j < 8*ElemOrder; j++) {
    //  //    Kcapa[N/2-4*ElemOrder+i][N/2-4*ElemOrder+j] = Rcapa_inv[N/2-4*ElemOrder+i][N/2-4*ElemOrder+j];
    //  //    Kelas[N/2-4*ElemOrder+i][N/2-4*ElemOrder+j] = Relas_inv[N/2-4*ElemOrder+i][N/2-4*ElemOrder+j];
    //  //  }
    //  //}

      Rcapa = 0;
      Relas = 0;
      for (Long i = 0; i < Rcapa.Dim(0); i++) Rcapa[i][i] = 1;
      for (Long i = 0; i < Relas.Dim(0); i++) Relas[i][i] = 1;
    }
  }
  convert_mat<double>(Kcapa).Write("Kcapa.mat");
  convert_mat<double>(Kelas).Write("Kelas.mat");

  if (0) { ///////////////
    Matrix<Real> P(N,N);
    P = 0;
    for (Long i = 0; i < N/2-4*ElemOrder; i++) P[i][i] = 1;
    for (Long i = 0; i < 8*ElemOrder; i++) P[N-8*ElemOrder+i][N/2-4*ElemOrder+i] = 1;
    for (Long i = 0; i < N/2-4*ElemOrder; i++) P[N/2-4*ElemOrder+i][N/2+4*ElemOrder+i] = 1;

    //(P*Kcapa*P.Transpose()).Write("Kcapa.mat");
    //(P*Kelas*P.Transpose()).Write("Kelas.mat");
    //Kcapa.Write("Kcapa.mat");
    //Kelas.Write("Kelas.mat");
    //Rcapa.Write("Rcapa.mat");
    //Relas.Write("Relas.mat");
  }

  auto SolveCapacitance = [&Kcapa, &Rcapa](const Vector<DiscType>& disc_lst, const Real tol, const Long max_iter) {
    Vector<Real> U0; // Dirichlet boundary conditions
    for (Long i = 0; i < disc_lst.Dim(); i++) { // Set U0
      const Long N = disc_lst[i].NodeCount();
      for (Long j = 0; j < N; j++) {
        U0.PushBack(i+1); // U0 = i for i-th disc
      }
    }

    // Setup boundary integral operator
    auto BIOp = [&Kcapa,&Rcapa](Vector<Real>* Ax, const Vector<Real>& x) {
      const Long N = x.Dim();
      if (Ax->Dim() != N) Ax->ReInit(N);
      Matrix<Real> U(1, N, Ax->begin(), false);
      Matrix<Real>::GEMM(U, Matrix<Real>(1, N, (Iterator<Real>)x.begin(), false), Kcapa);
    };

    Vector<Real> sigma, sigma_;
    if (max_iter == 0) {
      Matrix<Real> Kcapa_inv = Matrix<Real>(Kcapa).pinv(tol);
      Matrix<Real> sigma__ = Matrix<Real>(1,U0.Dim(),U0.begin(),false) * Kcapa_inv;
      sigma_.ReInit(sigma__.Dim(1), sigma__.begin());
    } else {
      ParallelSolver<Real> solver; // Linear solver (GMRES)
      solver(&sigma_, BIOp, U0, tol, max_iter);
    }
    { // Set sigma <-- sigma * R
      const Long N = sigma_.Dim();
      sigma.ReInit(N);
      Matrix<Real> Msigma(1,N, sigma.begin(), false);
      Msigma = Matrix<Real>(1,N, sigma_.begin(), false) * Rcapa;
    }

    // Return charge on each disk
    Vector<Real> Q;
    Long offset = 0;
    for (auto& disc : disc_lst) {
      const Long N = disc.NodeCount();
      Vector<Real> sigma_(N, sigma.begin() + offset, false);

      Vector<Real> q;
      disc.BoundaryIntegralDirect(q, sigma_);
      Q.PushBack(q[0]);

      offset += N;
    }
    return Q;
  };

  auto SolveElastance = [&Kelas, &Relas](const Vector<DiscType>& disc_lst, const Vector<Real> Q, const Real tol, const Long max_iter) {
    using KerSL = Laplace2D_FxU;
    using KerDL = Laplace2D_DxU;

    Vector<Real> Xt;
    GetGeom(disc_lst, &Xt);

    Vector<Real> nu, S_nu;
    SCTL_ASSERT(disc_lst.Dim() == Q.Dim());
    for (Long i = 0; i < disc_lst.Dim(); i++) {
      Real nu0 = Q[i] / (2*const_pi<Real>()*disc_lst[i].Radius());
      for (Long j = 0; j < disc_lst[i].NodeCount(); j++) {
        nu.PushBack(nu0);
      }
    }
    LayerPotential<KerSL>(S_nu, disc_lst, Xt, nu, tol); // Snu <-- S[nu]

    // Setup boundary integral operator
    auto BIOp = [&Kelas,&Relas](Vector<Real>* Ax, const Vector<Real>& x) {
      const Long N = x.Dim();
      if (Ax->Dim() != N) Ax->ReInit(N);
      Matrix<Real> U(1, N, Ax->begin(), false);
      Matrix<Real>::GEMM(U, Matrix<Real>(1, N, (Iterator<Real>)x.begin(), false), Kelas);
    };

    Vector<Real> sigma, sigma_; // unknown density
    ParallelSolver<Real> solver; // Linear solver (GMRES)
    if (max_iter == 0) {
      Matrix<Real> Kelas_inv = Matrix<Real>(Kelas).pinv(tol);
      Matrix<Real> sigma__ = Matrix<Real>(1,S_nu.Dim(),S_nu.begin(),false) * Kelas_inv * (-1);
      sigma_.ReInit(sigma__.Dim(1), sigma__.begin());
    } else {
      solver(&sigma_, BIOp, S_nu*(-1), tol, max_iter);
    }
    { // Set sigma <-- sigma * R
      const Long N = sigma_.Dim();
      sigma.ReInit(N);
      Matrix<Real> Msigma(1,N, sigma.begin(), false);
      Msigma = Matrix<Real>(1,N, sigma_.begin(), false) * Relas;
    }

    Vector<Real> D_sigma;
    LayerPotential<KerDL>(D_sigma, disc_lst, Xt, sigma, tol); // D_sigma <-- D[sigma]
    const Vector<Real> U = 0.5*sigma + D_sigma + S_nu;

    // print error
    Real max_err = 0;
    max_err = std::max<Real>(max_err, fabs(U[0]-1));
    max_err = std::max<Real>(max_err, fabs(U[U.Dim()-1]-2));
    return max_err;
  };

  sctl::Profile::Tic("Capacitance");
  const auto Q = SolveCapacitance(disc_lst, tol, max_iter);
  sctl::Profile::Toc();

  if (true || eps == atoreal<Real>("1e-10")) { // Print error ///////////////////////////////////
    Vector<Real> Q0;
    Q0.PushBack(atoreal<Real>("-2.720986561362458465987221035873359e+5"));
    Q0.PushBack(atoreal<Real>(" 2.720411531370425075198331580246059e+5"));

    const auto Qerr = Q/Q0-1;
    std::cout<<Qerr[0]<<' '<<Qerr[1]<<'\n';
    std::cout<<std::setprecision(34)<<Q[0]<<' '<<Q[1]<<'\n';
  }

  sctl::Profile::Tic("Elastance");
  Real err = SolveElastance(disc_lst, Q, tol, max_iter);
  sctl::Profile::Toc();

  return err;
}






template <class Real> Matrix<Real> BlockMatV(const Matrix<Real>& M0, const Matrix<Real>& M1) {
  const Long N0[2] = {M0.Dim(0), M1.Dim(0)};
  const Long N1 = {std::max<Long>(M0.Dim(1), M1.Dim(1))};
  Matrix<Real> M(N0[0]+N0[1], N1); M = 0;

  for (Long i = 0; i < M0.Dim(0); i++) {
    for (Long j = 0; j < M0.Dim(1); j++) {
      M[i][j] = M0[i][j];
    }
  }
  for (Long i = 0; i < M1.Dim(0); i++) {
    for (Long j = 0; j < M1.Dim(1); j++) {
      M[N0[0]+i][j] = M1[i][j];
    }
  }
  return M;
}
template <class Real> Matrix<Real> BlockMatH(const Matrix<Real>& M0, const Matrix<Real>& M1) {
  const Long N0 = {std::max<Long>(M0.Dim(0), M1.Dim(0))};
  const Long N1[2] = {M0.Dim(1), M1.Dim(1)};
  Matrix<Real> M(N0, N1[0]+N1[1]); M = 0;

  for (Long i = 0; i < M0.Dim(0); i++) {
    for (Long j = 0; j < M0.Dim(1); j++) {
      M[i][j] = M0[i][j];
    }
  }
  for (Long i = 0; i < M1.Dim(0); i++) {
    for (Long j = 0; j < M1.Dim(1); j++) {
      M[i][N1[0]+j] = M1[i][j];
    }
  }
  return M;
}
template <class Real> Matrix<Real> BlockMat(const Matrix<Real>& M00, const Matrix<Real>& M01, const Matrix<Real>& M10, const Matrix<Real>& M11) {
  const Long N0[2] = {std::max<Long>(M00.Dim(0), M01.Dim(0)),
                      std::max<Long>(M10.Dim(0), M11.Dim(0))};
  const Long N1[2] = {std::max<Long>(M00.Dim(1), M10.Dim(1)),
                      std::max<Long>(M01.Dim(1), M11.Dim(1))};
  Matrix<Real> M(N0[0]+N0[1], N1[0]+N1[1]); M = 0;

  for (Long i = 0; i < M00.Dim(0); i++) {
    for (Long j = 0; j < M00.Dim(1); j++) {
      M[i][j] = M00[i][j];
    }
  }
  for (Long i = 0; i < M01.Dim(0); i++) {
    for (Long j = 0; j < M01.Dim(1); j++) {
      M[i][N1[0]+j] = M01[i][j];
    }
  }
  for (Long i = 0; i < M10.Dim(0); i++) {
    for (Long j = 0; j < M10.Dim(1); j++) {
      M[N0[0]+i][j] = M10[i][j];
    }
  }
  for (Long i = 0; i < M11.Dim(0); i++) {
    for (Long j = 0; j < M11.Dim(1); j++) {
      M[N0[0]+i][N1[0]+j] = M11[i][j];
    }
  }
  return M;
}

template <class Real> struct DiscList {
  Vector<Real> R;
  Vector<Real> X;
  Vector<Vector<Real>> theta0, theta1;
};
template <class Real> DiscList<Real> GetAdapDisc(const Real R, const Real eps, const Long adap_depth, const Long disc_idx = -1) {
  auto theta_adap0 = [](Long depth) {
    Vector<Real> theta;
    for (Long i = 0; i < 20; i++) theta.PushBack(1+i*0.5);
    for (Long i = 0; i <=depth; i++) theta.PushBack(12-pow<Real>(0.5, i));
    theta.PushBack(12);
    for (Long i = 0; i <=depth; i++) theta.PushBack(12+pow<Real>(0.5, depth-i));
    return theta * const_pi<Real>()/6;
  }(adap_depth);
  auto theta_adap1 = [](Long depth) {
    Vector<Real> theta;
    for (Long i = 0; i <=depth; i++) theta.PushBack(6-pow<Real>(0.5, i));
    theta.PushBack(6);
    for (Long i = 0; i <=depth; i++) theta.PushBack(6+pow<Real>(0.5, depth-i));
    for (Long i = 0; i < 20; i++) theta.PushBack(7+(i+1)*0.5);
    return theta * const_pi<Real>()/6;
  }(adap_depth);

  DiscList<Real> disc_lst;
  if (disc_idx == 0 || disc_idx == -1) { // disc0
    disc_lst.R.PushBack(R);
    disc_lst.X.PushBack(-R-eps/2);
    disc_lst.X.PushBack( (Real)0);
    disc_lst.theta0.PushBack(Vector<Real>(theta_adap0.Dim()-1, theta_adap0.begin()+0));
    disc_lst.theta1.PushBack(Vector<Real>(theta_adap0.Dim()-1, theta_adap0.begin()+1));
  }
  if (disc_idx == 1 || disc_idx == -1) { // disc1
    disc_lst.R.PushBack(R);
    disc_lst.X.PushBack( R+eps/2);
    disc_lst.X.PushBack( (Real)0);
    disc_lst.theta0.PushBack(Vector<Real>(theta_adap1.Dim()-1, theta_adap1.begin()+0));
    disc_lst.theta1.PushBack(Vector<Real>(theta_adap1.Dim()-1, theta_adap1.begin()+1));
  }
  return disc_lst;
}
template <class Real> void CheckSurf(Vector<Real>* X0, Vector<Real>* X1, const Real R, Real eps) {
  const Vector<Real> nds0 = LegQuadRule<Real>::ComputeNds(std::max<Long>(1,ElemOrder));

  Vector<Real> theta0, theta1;
  for (Long i = 0; i <= 20; i++) theta0.PushBack(1+i*0.5);
  for (Long i = 0; i <= 20; i++) theta1.PushBack(7+i*0.5);
  theta0 = theta0 * const_pi<Real>()/6;
  theta1 = theta1 * const_pi<Real>()/6;

  Vector<Real> X0_, X1_;
  const Real tol = sqrt<Real>(machine_eps<Real>());
  for (Long i = 0; i < theta0.Dim()-1; i++) {
    for (Long j = 0; j < ElemOrder; j++) {
      { // Set X0
        const Real theta = theta0[i] + (theta0[i+1]-theta0[i]) * nds0[j];
        const Real x = -R-eps/2 + R*cos<Real>(theta) * (1+tol);
        const Real y =            R*sin<Real>(theta) * (1+tol);
        X0_.PushBack(x);
        X0_.PushBack(y);
      }
      { // Set X0
        const Real theta = theta0[i] + (theta0[i+1]-theta0[i]) * nds0[j];
        const Real x = -R-eps/2 + R*cos<Real>(theta) * (1-tol);
        const Real y =            R*sin<Real>(theta) * (1-tol);
        X0_.PushBack(x);
        X0_.PushBack(y);
      }
      { // Set X1
        const Real theta = theta1[i] + (theta1[i+1]-theta1[i]) * nds0[j];
        const Real x = R+eps/2 + R*cos<Real>(theta) * (1+tol);
        const Real y =           R*sin<Real>(theta) * (1+tol);
        X1_.PushBack(x);
        X1_.PushBack(y);
      }
      { // Set X1
        const Real theta = theta1[i] + (theta1[i+1]-theta1[i]) * nds0[j];
        const Real x = R+eps/2 + R*cos<Real>(theta) * (1-tol);
        const Real y =           R*sin<Real>(theta) * (1-tol);
        X1_.PushBack(x);
        X1_.PushBack(y);
      }
    }
  }

  const Long Nc = 1000;
  { // Set check surface
    const Real pi_6 = const_pi<Real>()/6;
    const Real R_sin_pi_6 = R*sin<Real>(pi_6);
    const Real R_cos_pi_6 = R*cos<Real>(pi_6);
    const Real RR = sqrt<Real>(R_sin_pi_6 * R_sin_pi_6 + (R+eps/2-R_cos_pi_6) * (R+eps/2-R_cos_pi_6));
    for (Long i = 0; i < Nc; i++) {
      const Real theta = i*2*const_pi<Real>()/Nc;
      const Real x = RR * cos<Real>(theta);
      const Real y = RR * sin<Real>(theta);
      X0_.PushBack(x);
      X0_.PushBack(y);
      X1_.PushBack(x);
      X1_.PushBack(y);
    }
  }

  if (X0 != nullptr) (*X0) = X0_;
  if (X1 != nullptr) (*X1) = X1_;
}
template <class Real> void CheckSurf_(Vector<Real>* X0, Vector<Real>* X1, const Real R, Real eps) {
  const Vector<Real> nds0 = LegQuadRule<Real>::ComputeNds(std::max<Long>(1,ElemOrder));

  //Vector<Real> theta0, theta1;
  //for (Long i = 0; i <= 20; i++) theta0.PushBack(1+i*0.5);
  //for (Long i = 0; i <= 20; i++) theta1.PushBack(7+i*0.5);
  //theta0 = theta0 * const_pi<Real>()/6;
  //theta1 = theta1 * const_pi<Real>()/6;

  Vector<Real> X0_, X1_;
  //const Real tol = sqrt<Real>(machine_eps<Real>());
  //for (Long i = 0; i < theta0.Dim()-1; i++) {
  //  for (Long j = 0; j < ElemOrder; j++) {
  //    { // Set X0
  //      const Real theta = theta0[i] + (theta0[i+1]-theta0[i]) * nds0[j];
  //      const Real x = -R-eps/2 + R*cos<Real>(theta) * (1+tol);
  //      const Real y =            R*sin<Real>(theta) * (1+tol);
  //      X0_.PushBack(x);
  //      X0_.PushBack(y);
  //    }
  //    { // Set X0
  //      const Real theta = theta0[i] + (theta0[i+1]-theta0[i]) * nds0[j];
  //      const Real x = -R-eps/2 + R*cos<Real>(theta) * (1-tol);
  //      const Real y =            R*sin<Real>(theta) * (1-tol);
  //      X0_.PushBack(x);
  //      X0_.PushBack(y);
  //    }
  //    { // Set X1
  //      const Real theta = theta1[i] + (theta1[i+1]-theta1[i]) * nds0[j];
  //      const Real x = R+eps/2 + R*cos<Real>(theta) * (1+tol);
  //      const Real y =           R*sin<Real>(theta) * (1+tol);
  //      X1_.PushBack(x);
  //      X1_.PushBack(y);
  //    }
  //    { // Set X1
  //      const Real theta = theta1[i] + (theta1[i+1]-theta1[i]) * nds0[j];
  //      const Real x = R+eps/2 + R*cos<Real>(theta) * (1-tol);
  //      const Real y =           R*sin<Real>(theta) * (1-tol);
  //      X1_.PushBack(x);
  //      X1_.PushBack(y);
  //    }
  //  }
  //}

  const Long Nc = 1000;
  { // Set check surface
    const Real pi_6 = const_pi<Real>()/6;
    const Real R_sin_pi_6 = R*sin<Real>(pi_6);
    const Real R_cos_pi_6 = R*cos<Real>(pi_6);
    const Real RR = sqrt<Real>(R_sin_pi_6 * R_sin_pi_6 + (R+eps/2-R_cos_pi_6) * (R+eps/2-R_cos_pi_6));
    for (Long i = 0; i < Nc; i++) {
      const Real theta = i*2*const_pi<Real>()/Nc;
      const Real x = RR * cos<Real>(theta);
      const Real y = RR * sin<Real>(theta);
      X0_.PushBack(x);
      X0_.PushBack(y);
      X1_.PushBack(x);
      X1_.PushBack(y);
    }
  }

  if (X0 != nullptr) (*X0) = X0_;
  if (X1 != nullptr) (*X1) = X1_;
}


template <class Real> Matrix<Real> FourierMatrix(const Long Nfourier, const Vector<Real>& theta0, const Vector<Real>& theta1, const Long Order) {
  const Vector<Real> nds0 = LegQuadRule<Real>::ComputeNds(std::max<Long>(1,Order));
  Vector<Real> theta(theta0.Dim()*Order);
  if (Order == 0) {
    theta = theta0;
  } else {
    for (Long i = 0; i < theta0.Dim(); i++) {
      for (Long j = 0; j < Order; j++) {
        theta[i*Order+j] = theta0[i] + (theta1[i]-theta0[i]) * nds0[j];
      }
    }
  }
  Matrix<Real> M(Nfourier, theta.Dim());
  for (Long i = 0; i < Nfourier; i++) {
    const Real scal = (i==Nfourier-1 && Nfourier%2==0 ? 1 : sqrt<Real>((Real)2));
    for (Long j = 0; j < theta.Dim(); j++) {
      if (i == 0) {
        M[i][j] = 1;
      } else if (i % 2 == 1) {
        const Long k = (i+1)/2;
        M[i][j] = cos<Real>(k*theta[j]) * scal;
      } else {
        const Long k = i/2;
        M[i][j] = sin<Real>(k*theta[j]) * scal;
      }
    }
  }
  return M;
}
template <class Real> Matrix<Real> FourierMatrix(const Long Nfourier, const Long Npanel, const Long Order) {
  Vector<Real> theta0(Npanel);
  Vector<Real> theta1(Npanel);
  for (Long i = 0; i < Npanel; i++) {
    theta0[i] = (i+0)*2*const_pi<Real>()/Npanel;
    theta1[i] = (i+1)*2*const_pi<Real>()/Npanel;
  }
  const Matrix<Real> M = FourierMatrix(Nfourier, theta0, theta1, Order);
  if (0) if (Nfourier == Npanel && Order == 0) { // Check that the matrix is orthonormal
    Real err = 0;
    Matrix<Real> I0 = M.Transpose() * M * (1/(Real)Nfourier);
    Matrix<Real> I1 = M * M.Transpose() * (1/(Real)Nfourier);
    for (Long i = 0; i < Nfourier; i++) {
      I0[i][i] -= 1;
      I1[i][i] -= 1;
    }
    for (Long i = 0; i < Nfourier*Nfourier; i++) {
      err = std::max<Real>(err, fabs(I0[0][i]));
      err = std::max<Real>(err, fabs(I1[0][i]));
    }
    SCTL_ASSERT(err < 4 * Nfourier * machine_eps<Real>());
  }
  return M;
}

template <class Real> Matrix<Real> Fourier2Nds(const Long N, const DiscList<Real> disc_lst) {
  Matrix<Real> Mf0 = FourierMatrix<Real>(N, N, 0);
  Matrix<Real> Mf0_inv = Mf0.Transpose() * (1/(Real)N); // Matrix<Real>(Mf0).pinv(machine_eps<Real>());

  Matrix<Real> M;
  const Long N_disc = disc_lst.R.Dim();
  for (Long i = 0; i < N_disc; i++) {
    Matrix<Real> Mf1 = FourierMatrix<Real>(N, disc_lst.theta0[i], disc_lst.theta1[i], ElemOrder);
    M = BlockMat(M, Matrix<Real>(), Matrix<Real>(), Mf0_inv * Mf1);
  }
  return M;
}
template <class Real> Matrix<Real> Nds2Fourier(DiscList<Real> disc_lst, const Long N) {
  static const Vector<Real> nds0 = LegQuadRule<Real>::ComputeNds(ElemOrder);
  const Long N_disc = disc_lst.R.Dim();

  Long Nnds = 0;
  for (const auto& theta0 : disc_lst.theta0) Nnds += theta0.Dim();
  Matrix<Real> M(Nnds*ElemOrder, N_disc*N);
  M.SetZero();

  Long nds_offset = 0;
  for (Long i = 0; i < N_disc; i++) {
    Vector<Real>& theta0 = disc_lst.theta0[i];
    Vector<Real>& theta1 = disc_lst.theta1[i];
    const Long Npanels = theta0.Dim();
    SCTL_ASSERT(theta1.Dim() == Npanels);
    for (Long j = 0; j < Npanels; j++) {
      while (theta1[j] < 0) {
        theta0[j] += 2*const_pi<Real>();
        theta1[j] += 2*const_pi<Real>();
      }
      while (theta1[j] > 2*const_pi<Real>()) {
        theta0[j] -= 2*const_pi<Real>();
        theta1[j] -= 2*const_pi<Real>();
      }
    }

    for (Long j = 0; j < N; j++) {
      Long idx = 0;
      Real theta = j*2*const_pi<Real>()/N;
      for (Long k = 0; k < Npanels; k++) {
        const Real theta_ = -2*const_pi<Real>();
        if (theta0[k] <= theta_ && theta_ <= theta1[k]) {
          theta = theta_;
          idx = k;
          break;
        }
        if (theta0[k] <= theta && theta <= theta1[k]) {
          idx = k;
          break;
        }
      }
      SCTL_ASSERT(theta0[idx] <= theta && theta <= theta1[idx]);

      StaticArray<Real,ElemOrder> interp_wts_;
      Vector<Real> interp_wts(ElemOrder, interp_wts_, false);
      StaticArray<Real,1> trg_node;
      trg_node[0] = {(theta-theta0[idx])/(theta1[idx]-theta0[idx])};
      LagrangeInterp<Real>::Interpolate(interp_wts, nds0, Vector<Real>(1, trg_node, false));
      for (Long k = 0; k < ElemOrder; k++) M[(nds_offset+idx)*ElemOrder+k][i*N+j] = interp_wts[k];
    }
    nds_offset += Npanels;
  }
  return M;
}
template <class Real> Matrix<Real> Nds2FourierProj(const DiscList<Real> disc_lst, const Long N) {
  Matrix<Real> Mf0 = FourierMatrix<Real>(N, N, 0);

  Matrix<Real> M;
  const Long N_disc = disc_lst.R.Dim();
  for (Long i = 0; i < N_disc; i++) {
    Matrix<Real> Mf1 = FourierMatrix<Real>(N, disc_lst.theta0[i], disc_lst.theta1[i], ElemOrder);
    Matrix<Real> Mf1_inv; // = Matrix<Real>(Mf1).pinv(machine_eps<Real>());
    { // set Mf1_inv = diag(W) * Mf1.Transpose() * scal // assumes Mf1 has orthogonal columns (but not orthonormal) wrt GL weights
      Vector<Real> W;
      using DiscType = Disc<Real,ElemOrder>;
      const DiscType disc = DiscType((Real)0, (Real)0, (Real)1, disc_lst.theta0[i], disc_lst.theta1[i]);
      disc.BoundaryIntegralWts(W);

      Mf1_inv = Mf1.Transpose();
      for (Long j = 0; j < Mf1_inv.Dim(0); j++) {
        for (Long k = 0; k < Mf1_inv.Dim(1); k++) {
          Mf1_inv[j][k] *= W[j]/(2*const_pi<Real>());
        }
      }
      if (N%2 == 0) {
        const Long  k = Mf1_inv.Dim(1)-1;
        for (Long j = 0; j < Mf1_inv.Dim(0); j++) Mf1_inv[j][k] *= 2;
      }
    }
    M = BlockMat(M, Matrix<Real>(), Matrix<Real>(), Mf1_inv * Mf0);
  }
  return M;
}


template <class Real> void CapaElasMatAdap(Matrix<Real>& Kcapa, Matrix<Real>& Kelas, const DiscList<Real>& disc_lst, const Real tol, const bool use_sqrtscal) {
  using DiscType = Disc<Real,ElemOrder>;

  const Long Ndisc = disc_lst.R.Dim();
  Vector<DiscType> disc_lst_(Ndisc);
  for (Long i = 0; i < Ndisc; i++) {
    disc_lst_[i] = DiscType(disc_lst.X[2*i+0], disc_lst.X[2*i+1], disc_lst.R[i], disc_lst.theta0[i], disc_lst.theta1[i]);
  }
  CapaElasOpMatrix(Kcapa, Kelas, disc_lst_, tol, use_sqrtscal);
}



template <class Ker, class Real> Matrix<Real> KerMatAdap(const DiscList<Real>& disc_lst, const Vector<Real>& Xt, const Real tol) {
  using DiscType = Disc<Real,ElemOrder>;

  const Long Ndisc = disc_lst.R.Dim();
  Vector<DiscType> disc_lst_(Ndisc);
  for (Long i = 0; i < Ndisc; i++) {
    disc_lst_[i] = DiscType(disc_lst.X[2*i+0], disc_lst.X[2*i+1], disc_lst.R[i], disc_lst.theta0[i], disc_lst.theta1[i]);
  }

  Vector<Real> Xt_(Xt.Dim(), (Iterator<Real>)Xt.begin(), false);
  if (!Xt_.Dim()) GetGeom(disc_lst_, &Xt_);

  Matrix<Real> M;
  LayerPotentialMatrix<Ker>(M, disc_lst_, Xt_, tol);
  return M;
}

template <class Real> void GetFourierEquivMat(Matrix<Real>& Mcapa, Matrix<Real>& Melas, const Matrix<Real>& Rcapa_, const Matrix<Real>& Relas_, const Long N, const Real R, const Real eps, const Long adap_depth, const Real tol,     const Matrix<Real>& Kcapa0) {
  const Matrix<Real> P1 = [&N]{
    const Long M = N/12;
    Matrix<Real> P(2*N, 2*(N-(2*M+1)));
    P.SetZero();
    for (Long i = 0; i < N-(2*M+1); i++) {
      P[M+1 + i][i] = 1;
    }
    for (Long i = 0; i < N/2-M-1; i++) {
      P[2*N-N/2+M+1 + i][N-(2*M+1) + i] = 1;
    }
    for (Long i = 0; i < N/2-M; i++) {
      P[N + i][2*(N-(2*M+1))-N/2+M + i] = 1;
    }
    return P;
  }();
  const Matrix<Real> P2 = [&N]{
    const Long M = N/12;
    Matrix<Real> P(2*N, 2*(2*M+1));
    P.SetZero();
    for (Long i = 0; i < M; i++) {
      P[N-M + i][i] = 1;
    }
    for (Long i = 0; i < M+1; i++) {
      P[i][M + i] = 1;
    }
    for (Long i = 0; i < 2*M+1; i++) {
      P[N+N/2-M + i][2*M+1 + i] = 1;
    }
    return P;
  }();
  const Matrix<Real> P = BlockMatH(P1,P2);
  P.template Write<double>((std::string(data_path) + "P.mat").c_str());

  Vector<Real> X0, X1;
  CheckSurf(&X0, &X1, R, eps);

  const auto disc_lst = GetAdapDisc(R, eps, adap_depth);
  const auto disc_lst0 = GetAdapDisc(R, eps, adap_depth, 0);
  const auto disc_lst1 = GetAdapDisc(R, eps, adap_depth, 1);

  const Matrix<Real> Mf2n = Fourier2Nds(N, disc_lst);
  const Matrix<Real> Mn2f = Nds2Fourier(disc_lst, N);

  Profile::Tic(std::to_string(__LINE__).c_str());
  const Matrix<Real> S0 = KerMatAdap<Laplace2D_FxU>(disc_lst0, X0, tol*1e-1);
  const Matrix<Real> S1 = KerMatAdap<Laplace2D_FxU>(disc_lst1, X1, tol*1e-1);
  const Matrix<Real> S = BlockMat(S0, Matrix<Real>(), Matrix<Real>(), S1);
  Profile::Toc();

  Profile::Tic(std::to_string(__LINE__).c_str());
  const Matrix<Real> D0 = KerMatAdap<Laplace2D_DxU>(disc_lst0, X0, tol*1e-1);
  const Matrix<Real> D1 = KerMatAdap<Laplace2D_DxU>(disc_lst1, X1, tol*1e-1);
  const Matrix<Real> D = BlockMat(D0, Matrix<Real>(), Matrix<Real>(), D1);
  Profile::Toc();

  Profile::Tic(std::to_string(__LINE__).c_str());
  {
    const auto Rcapa = Rcapa_*Mn2f + apply_inv((Rcapa_-Rcapa_*Mn2f*Mf2n) * (S+D), P2.Transpose()*Mf2n*(S+D), tol)* P2.Transpose(); // apply_inv(Rcapa_*(S+D), Mf2n*(S+D), tol);
    //const auto Rcapa = Rcapa_*Mn2f*P1*P1.Transpose() + apply_inv((Rcapa_-Rcapa_*Mn2f*P1*P1.Transpose()*Mf2n) * (S+D), P2.Transpose()*Mf2n*(S+D), tol)* P2.Transpose();
    //const auto Q = [P2,Mf2n,S,D]() {
    //  Matrix<Real> U, SS, Vt;
    //  (P2.Transpose()*Mf2n*(S+D)).SVD(U,SS,Vt);
    //  Real SSmax = 0;
    //  Vector<Real> UU;
    //  for (Long i = 0; i < SS.Dim(0); i++) {
    //    SSmax = std::max<Real>(SSmax, fabs(SS[i][i]));
    //  }
    //  for (Long j = 0; j < U.Dim(1); j++) {
    //    if (SS[j][j] < 1e-4*SSmax) {
    //      for (Long i = 0; i < U.Dim(0); i++) {
    //        UU.PushBack(U[i][j]);
    //      }
    //    }
    //  }
    //  return Matrix<Real>(UU.Dim()/U.Dim(0), U.Dim(0), UU.begin());
    //}();
    //Q.template Write<double>("data2/Q.mat");
    //std::cout<<Q.Dim(0)<<' '<<Q.Dim(1)<<'\n';

    ////const auto II = Matrix<Real>(Kcapa0).pinv() * P1*P1.Transpose() * Kcapa0;
    ////const auto X = II * Rcapa * P2*Q.Transpose();

    //const auto PP = [P1,Kcapa0]() {
    //  auto M = P1.Transpose()*Kcapa0;
    //  Matrix<Real> U,S,Vt;
    //  M.SVD(U,S,Vt);
    //  return Vt.Transpose() * Vt;
    //}();

    //const auto X = PP * (Rcapa*P2) * (Q.Transpose()*Q);
    //X.template Write<double>("data2/X.mat");

    Mcapa = Rcapa;// - X*P2.Transpose();

    //((P1.Transpose()*Kcapa0)*(Rcapa*P2)).template Write<double>("data2/E0.mat");
    ((P1.Transpose()*Kcapa0)*(Mcapa*P2)).template Write<double>("data2/E.mat");
  }
  if (0) {
    const auto Mf = FourierMatrix<Real>(N/3,N,0);
    const auto QQ = BlockMat(Mf.Transpose()*Mf, Matrix<Real>(), Matrix<Real>(), Mf.Transpose()*Mf) * (1/(Real)Mf.Dim(1));
    auto identity = [](const Long N) {
      Matrix<Real> M(N,N); M = 0;
      for (Long i = 0; i < N; i++) M[i][i] = 1;
      return M;
    };
    const auto Proj0 = [&P1,&Kcapa0,&QQ,N]() {
      const auto FourierPOU = [](const Long Nfourier, const Long Npanel, const Long Order, const Real theta_c) {
        const auto pou = [theta_c](Real theta) {
          return 1 - exp<Real>(-36 * pow<24>(0.5*(theta-theta_c))) - exp<Real>(-36 * pow<24>(0.5*(theta-(theta_c+2*const_pi<Real>()))));
        };
        const Vector<Real> nds0 = LegQuadRule<Real>::ComputeNds(std::max<Long>(1,Order));

        Vector<Real> theta0(Npanel);
        Vector<Real> theta1(Npanel);
        for (Long i = 0; i < Npanel; i++) {
          theta0[i] = (i+0)*2*const_pi<Real>()/Npanel;
          theta1[i] = (i+1)*2*const_pi<Real>()/Npanel;
        }

        Vector<Real> theta(theta0.Dim()*Order);
        if (Order == 0) {
          theta = theta0;
        } else {
          for (Long i = 0; i < theta0.Dim(); i++) {
            for (Long j = 0; j < Order; j++) {
              theta[i*Order+j] = theta0[i] + (theta1[i]-theta0[i]) * nds0[j];
            }
          }
        }
        Matrix<Real> M(Nfourier, theta.Dim());
        for (Long i = 0; i < Nfourier; i++) {
          const Real scal = (i==Nfourier-1 && Nfourier%2==0 ? 1 : sqrt<Real>((Real)2));
          for (Long j = 0; j < theta.Dim(); j++) {
            if (i == 0) {
              M[i][j] = 1;
            } else if (i % 2 == 1) {
              const Long k = (i+1)/2;
              M[i][j] = cos<Real>(k*theta[j]) * scal;
            } else {
              const Long k = i/2;
              M[i][j] = sin<Real>(k*theta[j]) * scal;
            }
            M[i][j] *= pou(theta[j]);
          }
        }
        return M;
      };
      const auto F_pou = BlockMat(FourierPOU(N/3,N,0,0), Matrix<Real>(), Matrix<Real>(), FourierPOU(N/3,N,0,const_pi<Real>()));
      F_pou.template Write<double>("data2/F_pou.mat");

      Matrix<Real> U, S, Vt;
      (QQ*Kcapa0).SVD(U,S,Vt);
      //(P1.Transpose()*Kcapa0).SVD(U,S,Vt);
      //(F_pou*Kcapa0).SVD(U,S,Vt);
      Real Smax = 0;
      for (Long i = 0; i < S.Dim(0); i++) {
        Smax = std::max<Real>(Smax, fabs(S[i][i]));
      }
      for (Long i = 0; i < Vt.Dim(0); i++) {
        for (Long j = 0; j < Vt.Dim(1); j++) {
          if (S[i][i] < 1e-10*Smax) {
            Vt[i][j] = 0;
          }
        }
      }
      return Vt.Transpose() * Vt;
    }();
    const auto Proj1 = identity(Proj0.Dim(0)) - Proj0;
    Proj0.template Write<double>((std::string(data_path) + "Proj0.mat").c_str());

    const auto Rcapa0 = Matrix<Real>(Kcapa0).pinv(tol);
    const auto Rcapa = Rcapa_*Mn2f + apply_inv((Rcapa_-Rcapa_*Mn2f*Mf2n) * (S+D), P2.Transpose()*Mf2n*(S+D), tol)* P2.Transpose(); // apply_inv(Rcapa_*(S+D), Mf2n*(S+D), tol);

    const auto U0 = Proj0 * Rcapa_     *(S+D);
    const auto U1 = Proj0 * Rcapa0*Mf2n*(S+D);
    const auto U2 = Proj0 * Rcapa *Mf2n*(S+D);
    //const auto U0 = QQ*Kcapa0 * Rcapa_     *(S+D);
    //const auto U1 = QQ*Kcapa0 * Rcapa0*Mf2n*(S+D);
    //const auto U2 = QQ*Kcapa0 * Rcapa *Mf2n*(S+D);
    //U0.template Write("data2/U0.mat");
    //U1.template Write("data2/U1.mat");
    //U2.template Write("data2/U2.mat");
    const auto max_rel_err = [](Matrix<Real> U0, const Matrix<Real>& U) {
      const auto E = U0 - U;
      Real rel_err = 0;
      for (Long i = 0; i < E.Dim(0); i++) {
        Real max_err = 0, max_val = 0;
        for (Long j = 0; j < E.Dim(1); j++) {
          max_err = std::max<Real>(max_err, fabs(E[i][j]));
          max_val = std::max<Real>(max_val, fabs(U0[i][j]));
        }
        rel_err = std::max<Real>(rel_err, max_err/max_val);
      }
      return rel_err;
    };
    std::cout<<"Compression accuracy = "<<max_rel_err(U0, U1)<<'\n';
    std::cout<<"Compression accuracy = "<<max_rel_err(U0, U2)<<'\n';

    //Mcapa = Proj0*Rcapa0 + Proj1*Rcapa;  //*Proj1*apply_inv((Rcapa_-Mcapa*Mf2n) * (S+D), Mf2n*(S+D), tol);
    //Mcapa = Rcapa0;

    Mcapa = Rcapa; // Rcapa0
    Mcapa = Mcapa + QQ*(Rcapa-Mcapa);
    Mcapa = Mcapa + Proj0*(Rcapa0-Mcapa);

    Mcapa = Rcapa0;
    Mcapa = Mcapa + QQ*(Rcapa-Mcapa);
    Mcapa = Mcapa + (identity(600)-QQ)*Proj0*(Rcapa0-Mcapa);
    // Orthogonality of QQ and Proj is the issue

    Mcapa = Rcapa0;
    Profile::Tic(std::to_string(__LINE__).c_str());
    for (Long i = 0; i < 1; i++) {
      Mcapa = Mcapa + QQ*(Rcapa - Mcapa);
      Mcapa = Mcapa + Proj0*(Rcapa0 - Mcapa);
    }
    Profile::Toc();

    auto Z = QQ*(identity(600)-Proj0);
    { // Construct projection to span(Z)
      Matrix<Real> U,S,Vt;
      Z.SVD(U,S,Vt);
      Real Smax = 0;
      for (Long i = 0; i < S.Dim(0); i++) {
        Smax = std::max<Real>(Smax, fabs(S[i][i]));
      }
      for (Long i = 0; i < Vt.Dim(0); i++) {
        for (Long j = 0; j < Vt.Dim(1); j++) {
          if (S[i][i] < 1e-6*Smax) {
            Vt[i][j] = 0;
          }
        }
      }
      Z = Vt.Transpose()*Vt;
    }
    Z.template Write<double>("data2/Z.mat");
    QQ.template Write<double>("data2/QQ.mat");
    Proj0.template Write<double>("data2/Proj0.mat");
    Mcapa = Rcapa0;
    Mcapa = Mcapa + Z * (Rcapa - Mcapa);
  }

  //Mcapa = Rcapa_*Mn2f + apply_inv((Rcapa_-Rcapa_*Mn2f*Mf2n) * (S+D), P2.Transpose()*Mf2n*(S+D), tol)* P2.Transpose();
  //Melas = Relas_*Mn2f + apply_inv((Relas_-Relas_*Mn2f*Mf2n) *    D , P2.Transpose()*Mf2n*   D , tol)* P2.Transpose();

  //Mcapa = Rcapa_*Mn2f*P1*P1.Transpose() + apply_inv(Rcapa_*(S+D), Mf2n*(S+D), tol*0+1e-7)*P2*P2.Transpose();
  //Melas = Relas_*Mn2f*P1*P1.Transpose() + apply_inv(Relas_*   D , Mf2n*   D , tol*0+1e-7)*P2*P2.Transpose();

  //Mcapa = apply_inv(Rcapa_*(S+D), Mf2n*(S+D), tol*0+1e-7);
  //Melas = apply_inv(Relas_*   D , Mf2n*   D , tol*0+1e-7);
  Profile::Toc();

  // K(1:498,:) * (R*Mn2f)(:,499:end) = 0; // is this true already??
  // K(1:498,:) * M(:,499:end) = 0; // need to enforce
  {
    //const auto KK = P.Transpose()*Kcapa0*P;
    //const auto RR = P.Transpose()*Rcapa_*Mn2f*P;
    //const auto EE = KK*RR;
    //RR.template Write<double>((std::string(data_path) + "RR.mat").c_str());
    //KK.template Write<double>((std::string(data_path) + "KK.mat").c_str());
    //EE.template Write<double>((std::string(data_path) + "EE.mat").c_str());
    //Real err = 0;
    //for (Long i = 0; i < 498; i++) {
    //  for (Long j = 498; i < 600; i++) {
    //    err = std::max<Real>(err, fabs(EE[i][j]));
    //  }
    //}
    //std::cout<<err<<'\n';

    if (0) {
      Matrix<Real> Mf = FourierMatrix<Real>(N/3,N,0);
      Mf = BlockMat(Mf, Matrix<Real>(), Matrix<Real>(), Mf);

      Matrix<Real> U, S, Vt;
      (Mf*Kcapa0).SVD(U,S,Vt);
      auto Proj = Vt.Transpose() * Vt;
      //Mcapa = Proj * Rcapa_*Mn2f + (Mcapa-Proj*Mcapa);
      Mcapa = Proj * Matrix<Real>(Kcapa0).pinv(tol)  +  (Mcapa-Proj*Mcapa);

      // Kcapa0.pinv() = Vt.Transpose()*S.pinv()*U.Transpose()*Mf
    }

    if (0) {
      Matrix<Real> U, S, Vt;
      (P1.Transpose()*Kcapa0).SVD(U,S,Vt);
      auto Proj = Vt.Transpose() * Vt;
      //Mcapa = Proj * Rcapa_*Mn2f  +  (Mcapa-Proj*Mcapa);
      //Mcapa = Proj * Matrix<Real>(Kcapa0).pinv(tol)  +  (Mcapa-Proj*Mcapa);

      Mcapa = Proj * Mcapa  +  (Mcapa-Proj*Mcapa);
    }

    //Mcapa = Proj*Mcapa*P2*P2.Transpose() + Mcapa*P1*P1.Transpose();
  }



  Matrix<Real> sigma(1,600); sigma = 0;
  const Matrix<Real> Mf = FourierMatrix<Real>(N,N,0);
  for(Long i = 0; i < N; i++) {
    sigma[0][i] = Mf[1][i];
    //sigma[0][i] = 0;
    //for(Long j = 0; j < 600; j++) {
    //  sigma[0][i] += Kcapa0[j][i];
    //}
  }
  const auto u = sigma*Kcapa0;

  //const auto sigma_ = u*(Rcapa_-Rcapa_*Mn2f*Mf2n)*Mn2f;
  const auto sigma_ = u*Mcapa;
  //const auto sigma_ = u*Rcapa_*Mn2f;

  const auto err = sigma_-sigma;
  Real max_err = 0;
  for (const auto& x : err) max_err = std::max<Real>(max_err, fabs(x));
  std::cout<<max_err<<'\n';

}






int main(int argc, char** argv) {
  using Real = long double; //QuadReal;
  sctl::Profile::Enable(true);

  const auto write_precomp = [](const Matrix<Real>& M, const std::string& name, const Long idx) {
    std::string fname = std::string("precomp/") + name + std::to_string(idx) + ".mat";
    M.Write<QuadReal>(fname.c_str());
  };
  const auto read_precomp = [](const std::string& name, const Long idx) {
    std::string fname = std::string("precomp/") + name + std::to_string(idx) + ".mat";
    Matrix<Real> M;
    M.Read<QuadReal>(fname.c_str());
    return M;
  };

  const Real R = 0.75;
  const Real tol = machine_eps<Real>();
  const Real eps = atoreal<Real>("1e-10");
  const Long Nf = 300;
  const Long adap_depth = 15;

  std::cout<<"R = "<<R<<'\n';
  std::cout<<"Disc separation = "<<eps<<'\n';
  std::cout<<"Precision = "<<sizeof(Real)*8<<" bits\n";
  std::cout<<"ElemOrder = "<<ElemOrder<<'\n';
  std::cout<<"FourierOrder = "<<Nf<<'\n';
  std::cout<<"tol = "<<tol<<'\n';
  std::cout<<"data path = "<<data_path<<'\n';

  const auto disc_lst = GetAdapDisc(R, eps, adap_depth);
  const Matrix<Real> Mf2n = Fourier2Nds(Nf, disc_lst);
  const Matrix<Real> Mn2f = Nds2Fourier(disc_lst, Nf);
  if (0) { // Print accuracy of piecewise polynomials for representing Nf Fourier modes
    // ElemOrder=32, Nf = 150, Accuracy = 1e-15
    // ElemOrder=32, Nf = 200, Accuracy = 1e-11
    // ElemOrder=32, Nf = 250, Accuracy = 1e-08
    // ElemOrder=32, Nf = 300, Accuracy = 1e-06
    //
    // ElemOrder=48, Nf = 200, Accuracy = 1e-24
    // ElemOrder=48, Nf = 250, Accuracy = 1e-19
    // ElemOrder=48, Nf = 300, Accuracy = 1e-16
    //
    // ElemOrder=64, Nf = 300, Accuracy = 1e-28
    // ElemOrder=64, Nf = 400, Accuracy = 1e-21
    // ElemOrder=64, Nf = 600, Accuracy = 1e-10
    //
    // ElemOrder=96, Nf = 400, Accuracy = 1e-31
    // ElemOrder=96, Nf = 600, Accuracy = 1e-30
    // ElemOrder=96, Nf = 800, Accuracy = 1e-19

    Real err = 0;
    Matrix<Real> MM;
    const auto Mfourier10 = FourierMatrix<Real>(Nf, Nf*10, 0).pinv() * FourierMatrix<Real>(Nf, Nf, 0);
    for (Long i = 0; i < 2; i++) MM = BlockMat(MM, Matrix<Real>(), Matrix<Real>(), Mfourier10);
    auto Mf2f = Mf2n * Nds2Fourier(disc_lst, Nf*10) * MM;

    Mf2f.Write((std::string(data_path) + "Mf2f.mat").c_str());
    for (Long i = 0; i < Mf2f.Dim(0); i++) Mf2f[i][i] -= 1;
    for (Long i = 0; i < Mf2f.Dim(0)*Mf2f.Dim(1); i++) err = std::max<Real>(err, fabs(Mf2f[0][i]));
    std::cout<<"Elem discr err = "<<err<<'\n';
  }

  Profile::Tic("ComputeK");
  Matrix<Real> Kcapa_, Kelas_;
  Kcapa_ = read_precomp("Kcapa_", Mn2f.Dim(0));
  Kelas_ = read_precomp("Kelas_", Mn2f.Dim(0));
  if (Kcapa_.Dim(0)*Kcapa_.Dim(1) == 0 || Kelas_.Dim(0)*Kelas_.Dim(1) == 0) {
    CapaElasMatAdap(Kcapa_, Kelas_, disc_lst, tol, false);
    write_precomp(Kcapa_, "Kcapa_", Mn2f.Dim(0));
    write_precomp(Kelas_, "Kelas_", Mn2f.Dim(0));
  }
  Profile::Toc();

  Profile::Tic("ComputeR");
  Matrix<Real> Rcapa_, Relas_;
  Rcapa_ = read_precomp("Rcapa_", Mn2f.Dim(0));
  Relas_ = read_precomp("Relas_", Mn2f.Dim(0));
  if (Rcapa_.Dim(0)*Rcapa_.Dim(1) == 0 || Relas_.Dim(0)*Relas_.Dim(1) == 0) {
    Rcapa_ = Matrix<Real>(Kcapa_).pinv(tol);
    Relas_ = Matrix<Real>(Kelas_).pinv(tol);
    write_precomp(Rcapa_, "Rcapa_", Mn2f.Dim(0));
    write_precomp(Relas_, "Relas_", Mn2f.Dim(0));
  }
  Profile::Toc();

  ///////////////////////////////////////////////////////////////////////////
  Mf2n.Write<double>((std::string(data_path) + "Mf2n.mat").c_str());
  Mn2f.Write<double>((std::string(data_path) + "Mn2f.mat").c_str());
  //Kcapa_.Write<double>((std::string(data_path) + "Kcapa_.mat").c_str());
  //Rcapa_.Write<double>((std::string(data_path) + "Rcapa_.mat").c_str());
  ///////////////////////////////////////////////////////////////////////////

  Profile::Tic("ComputeEquiv");
  Matrix<Real> Rcapa; // = Mf2n * Rcapa_ * Mcapa;
  Matrix<Real> Relas; // = Mf2n * Relas_ * Melas;
  GetFourierEquivMat(Rcapa, Relas, Mf2n*Rcapa_,  Mf2n*Relas_, Nf, R, eps, adap_depth, tol,       Mf2n * Kcapa_ * Mn2f);
  Rcapa.Write<double>((std::string(data_path) + "Rcapa.mat").c_str());
  //Relas.Write<double>((std::string(data_path) + "Relas.mat").c_str());
  Profile::Toc();

  auto check_error = [eps] (const Vector<Real>& Q) { // Print error
    if (eps == atoreal<Real>("1e-10")) {
      Vector<Real> Q0;
      Q0.PushBack(atoreal<Real>("-2.720986561362458465987221035873359e+5"));
      Q0.PushBack(atoreal<Real>(" 2.720411531370425075198331580246059e+5"));

      const auto Qerr = Q/Q0-1;
      std::cout<<"Capacitance Error = "<<Qerr[0]<<' '<<Qerr[1]<<'\n';
    }
    //std::cout<<std::setprecision(34)<<Q[0]<<' '<<Q[1]<<"\n\n";
  };
  { // Solve capacitance (fine mesh)
    const Long N = Rcapa_.Dim(0);
    const Long Ndisc = disc_lst.R.Dim();

    Matrix<Real> U0(1,N);
    for (Long i = 0; i < N/2; i++) U0[0][i] = 1;
    for (Long i = N/2; i < N; i++) U0[0][i] = 2;
    const Matrix<Real> sigma = U0 * Rcapa_;

    using DiscType = Disc<Real,ElemOrder>;
    Vector<DiscType> disc_lst_(Ndisc);
    for (Long i = 0; i < Ndisc; i++) {
      disc_lst_[i] = DiscType(disc_lst.X[2*i+0], disc_lst.X[2*i+1], disc_lst.R[i], disc_lst.theta0[i], disc_lst.theta1[i]);
    }

    Vector<Real> Q;
    Long offset = 0;
    for (auto& disc : disc_lst_) {
      const Long N = disc.NodeCount();
      Vector<Real> sigma_(N, (Iterator<Real>)sigma.begin() + offset, false);

      Vector<Real> q;
      disc.BoundaryIntegralDirect(q, sigma_);
      Q.PushBack(q[0]);

      offset += N;
    }
    check_error(Q);
  }
  { // Solve capacitance
    const Long N = Rcapa.Dim(0);
    const Long Ndisc = disc_lst.R.Dim();

    Matrix<Real> U0(1,N);
    for (Long i = 0; i < N/2; i++) U0[0][i] = 1;
    for (Long i = N/2; i < N; i++) U0[0][i] = 2;
    const Matrix<Real> sigma = U0 * Rcapa;

    Vector<Real> Q(Ndisc); Q = 0;
    for (Long k = 0; k < Ndisc; k++) {
      for (Long i = 0; i < N/Ndisc; i++) {
        Q[k] += 2*const_pi<Real>()*R * sigma[0][k*N/Ndisc + i] / (N/Ndisc);
      }
    }
    check_error(Q);
  }

  // Uncompressed K and R
  (Mf2n * Kcapa_ * Mn2f).Write<double>((std::string(data_path) + "Kcapa0.mat").c_str());
  //(Mf2n * Kelas_ * Mn2f).Write<double>((std::string(data_path) + "Kelas0.mat").c_str());
  (Mf2n * Rcapa_ * Mn2f).Write<double>((std::string(data_path) + "Rcapa0.mat").c_str());
  //(Mf2n * Relas_ * Mn2f).Write<double>((std::string(data_path) + "Relas0.mat").c_str());

  Matrix<Real>(Rcapa).pinv(tol).Write<double>((std::string(data_path) + "Kcapa.mat").c_str());
  //Matrix<Real>(Relas).pinv(tol).Write<double>((std::string(data_path) + "Kelas.mat").c_str());

  { // Print accuracy of Rcapa
    const auto max_rel_err = [](Matrix<Real> U0, const Matrix<Real>& U) {
      const auto E = U0 - U;
      Real rel_err = 0;
      for (Long i = 0; i < E.Dim(0); i++) {
        Real max_err = 0, max_val = 0;
        for (Long j = 0; j < E.Dim(1); j++) {
          max_err = std::max<Real>(max_err, fabs(E[i][j]));
          max_val = std::max<Real>(max_val, fabs(U0[i][j]));
        }
        rel_err = std::max<Real>(rel_err, max_err/max_val);
      }
      return rel_err;
    };

    Profile::Tic("ComprAccuracy");
    Vector<Real> X, X0, X1;
    CheckSurf(&X0, &X1, R, eps);
    for (const auto& x : X0) X.PushBack(x);
    for (const auto& x : X1) X.PushBack(x);

    const auto S = KerMatAdap<Laplace2D_FxU>(disc_lst, X, tol);
    const auto D = KerMatAdap<Laplace2D_FxU>(disc_lst, X, tol);

    const auto Kcapa0 = (Mf2n * Kcapa_ * Mn2f);
    const auto Rcapa0 = (Mf2n * Rcapa_ * Mn2f);
    const auto Mf = FourierMatrix<Real>(75,300,0);
    const auto qq = BlockMat(Mf, Matrix<Real>(), Matrix<Real>(), Mf);

    //const auto U0 = qq * Mf2n*(S+D);
    //const auto U = qq*Mf2n*Kcapa_             * Rcapa_ *             (S+D); // 1e-27
    //const auto U = qq*Mf2n*Kcapa_             * Rcapa_ * Mn2f*Mf2n * (S+D); // 1e-27
    //const auto U = qq*Mf2n*Kcapa_ * Mn2f*Mf2n * Rcapa_             * (S+D); // 1e-11
    //const auto U = qq*Mf2n*Kcapa_ * Mn2f*Mf2n * Rcapa_ * Mn2f*Mf2n * (S+D); // 1e-11
    //const auto   U = qq*Mf2n*Kcapa_ * Mn2f * Matrix<Real>(Kcapa0).pinv(tol) * Mf2n*(S+D); // 1e-23
    //std::cout<<"Error ("<<__LINE__<<") = "<<max_rel_err(U0, U)<<'\n';

    {
      const auto U0 = qq * Mf2n*(S+D);
      const auto U = qq*Mf2n*Kcapa_*Mn2f * Rcapa * Mf2n*(S+D);
      std::cout<<"Compression accuracy = "<<max_rel_err(U0, U)<<'\n';
    }
    {
      const auto U0 = qq * Mf2n*Rcapa_ *      (S+D);
      const auto U  = qq *       Rcapa * Mf2n*(S+D);
      std::cout<<"Compression accuracy = "<<max_rel_err(U0, U)<<'\n';
    }

    Profile::Toc();
  }

  std::cout<<"cond(Rcapa) = "<<cond(Rcapa)<<'\n';
  //std::cout<<cond(Relas)<<'\n';

  Profile::print();

  return 0;
}

int main_(int argc, char** argv) {
  using Real = QuadReal;
  Comm::MPI_Init(&argc, &argv);
  sctl::Profile::Enable(true);

  const Real R = 0.75;

  const Real eps = atoreal<Real>("1e-10");
  const Real tol = 1e-20; //machine_eps<Real>();
  const Long max_iter = 0;

  std::cout<<"R = "<<R<<'\n';
  std::cout<<"Disc separation = "<<eps<<'\n';
  std::cout<<"Precision = "<<sizeof(Real)*8<<" bits\n";
  std::cout<<"ElemOrder = "<<ElemOrder<<'\n';
  std::cout<<"tol = "<<tol<<'\n';
  std::cout<<"Ngmres = "<<max_iter<<'\n';

  //while(true)
  {
    static Real max_err = 0, max_eps = 0;
    //const Real eps = pow<Real>(0.5, 30+6*drand48());
    const Real err = TestCompressPrecond(R, eps, tol, max_iter);
    if (err > max_err) {
      max_err = err;
      max_eps = eps;
    }
    std::cout<<"eps = "<<eps<<"    error = "<<err<<'\n';
    std::cout<<"eps = "<<max_eps<<"    error = "<<max_err<<'\n';
  }

  Comm::MPI_Finalize();
  return 0;

  const bool sqrt_scal = true;
  const bool unif_mesh = false;
  const bool use_gmres = false;
  const Long adap_depth = 17;

  Vector<Disc<Real,ElemOrder>> disc_lst(2);
  if (unif_mesh) { // Uniform mesh
    disc_lst[0] = Disc<Real,ElemOrder>(-R-eps/2, 0.0, R, 12);
    disc_lst[1] = Disc<Real,ElemOrder>( R+eps/2, 0.0, R, 12);
  } else { // Adaptive mesh
    auto theta_adap0 = [](Long depth) {
      Vector<Real> theta;
      for (Long i = 0; i < 10; i++) theta.PushBack(1+i);
      for (Long i = 0; i <=depth; i++) theta.PushBack(12-pow<Real>(0.5, i));
      theta.PushBack(12);
      for (Long i = 0; i <=depth; i++) theta.PushBack(12+pow<Real>(0.5, depth-i));
      return theta * const_pi<Real>()/6;
    }(adap_depth);
    auto theta_adap1 = [](Long depth) {
      Vector<Real> theta;
      for (Long i = 0; i <=depth; i++) theta.PushBack(6-pow<Real>(0.5, i));
      theta.PushBack(6);
      for (Long i = 0; i <=depth; i++) theta.PushBack(6+pow<Real>(0.5, depth-i));
      for (Long i = 0; i < 10; i++) theta.PushBack(7+(1+i));
      return theta * const_pi<Real>()/6;
    }(adap_depth);
    disc_lst[0] = Disc<Real,ElemOrder>(-R-eps/2, 0.0, R, Vector<Real>(theta_adap0.Dim()-1, theta_adap0.begin()+0), Vector<Real>(theta_adap0.Dim()-1, theta_adap0.begin()+1));
    disc_lst[1] = Disc<Real,ElemOrder>( R+eps/2, 0.0, R, Vector<Real>(theta_adap1.Dim()-1, theta_adap1.begin()+0), Vector<Real>(theta_adap1.Dim()-1, theta_adap1.begin()+1));
  }

  // Print parameters
  std::cout<<"R = "<<R<<'\n';
  std::cout<<"Disc separation = "<<eps<<'\n';
  std::cout<<"GMRES = "<<use_gmres<<'\n';
  std::cout<<"Sqrt-scaling = "<<sqrt_scal<<'\n';
  if (!unif_mesh) std::cout<<"Adap-depth = "<<adap_depth<<'\n';
  std::cout<<"Precision = "<<sizeof(Real)*8<<" bits\n";
  std::cout<<"ElemOrder = "<<ElemOrder<<'\n';
  std::cout<<"Nunknowns = "<<NodeCount(disc_lst)<<'\n';
  std::cout<<"tol = "<<tol<<'\n';
  std::cout<<"Ngmres = "<<max_iter<<'\n';

  sctl::Profile::Tic("Capacitance");
  const auto Q = SolveCapacitance<Real>(disc_lst, use_gmres, sqrt_scal, tol, max_iter);
  sctl::Profile::Toc();

  sctl::Profile::Tic("Elastance");
  SolveElastance(disc_lst, Q, use_gmres, sqrt_scal, tol, max_iter);
  sctl::Profile::Toc();

  sctl::Profile::print();
  return 0;
}

