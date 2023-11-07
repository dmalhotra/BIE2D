namespace sctl {

  template <class Real, Integer Order> DiscElastance<Real,Order>::DiscElastance(const Comm& comm_) : ICIP_Base(comm_), LaplaceSL_BIOp(LaplaceSL_Ker, false, this->comm), LaplaceSL_BIOp_near(LaplaceSL_Ker, false, this->comm), LaplaceSL_BIOp_far(LaplaceSL_Ker, false, this->comm), LaplaceDL_BIOp_near(LaplaceDL_Ker, false, this->comm), LaplaceDL_BIOp_far(LaplaceDL_Ker, false, this->comm), solver(this->comm, true) {}

  template <class Real, Integer Order> void DiscElastance<Real,Order>::Init(const Vector<Real>& Xc, const Real R, const Real tol, const ICIPType icip_type) {
    ICIP_Base::Init(Xc, R, tol, icip_type);

    LaplaceSL_BIOp     .AddElemList(this->disc_panels); // TODO: this can be just near2near
    LaplaceSL_BIOp_near.AddElemList(this->panels_near);
    LaplaceSL_BIOp_far .AddElemList(this->panels_far );
    LaplaceDL_BIOp_near.AddElemList(this->panels_near);
    LaplaceDL_BIOp_far .AddElemList(this->panels_far );

    LaplaceSL_BIOp_near.SetAccuracy(tol);
    LaplaceSL_BIOp_far .SetAccuracy(tol);
    LaplaceDL_BIOp_near.SetAccuracy(tol);
    LaplaceDL_BIOp_far .SetAccuracy(tol);

    LaplaceSL_BIOp_near.SetTargetCoord(this->Xfar);
    LaplaceSL_BIOp_far .SetTargetCoord(this->X   );
    LaplaceDL_BIOp_near.SetTargetCoord(this->Xfar);
    LaplaceDL_BIOp_far .SetTargetCoord(this->X   );
  }

  template <class Real, Integer Order> void DiscElastance<Real,Order>::Solve(Vector<Real>& V, const Vector<Real> Q, const Long gmres_max_iter) {
    const Long Ndisc = this->disc_panels.DiscCount();
    const Long N = this->disc_panels.Size() * Order;
    const Real R = this->disc_panels.DiscRadius();

    Vector<Real> nu; // boundary force
    { // Set nu
      nu.ReInit(N);
      Long offset = 0;
      for (Long i = 0; i < Ndisc; i++) {
        const Long N_ = this->disc_panels.SurfWts(i).Dim();
        const Real inv_circumf = 1/(2*const_pi<Real>()*R);
        Vector<Real> nu_(N_, nu.begin() + offset, false);
        nu_ = Q[i] * inv_circumf;
        offset += N_;
      }
    }

    Vector<Real> U0; // completion flow velocity
    LaplaceSL_BIOp.ComputePotential(U0, nu);

    Vector<Real> sigma, sigma_; // unknown density
    const auto LaplaceElastOp = [this](Vector<Real>* U, const Vector<Real> sigma) {
      this->ApplyBIOp(U, sigma);
    };
    solver(&sigma_, LaplaceElastOp, U0, this->tol_, gmres_max_iter);

    this->ApplyPrecond(&sigma, sigma_);

    if (V.Dim() != Ndisc) V.ReInit(Ndisc);
    { // get average of potential
      V = 0;
      Long offset = 0;
      for (Long j = 0; j < Ndisc; j++) {
        const auto& wts = this->disc_panels.SurfWts(j);
        const Long N_ = wts.Dim();
        for (Long i = 0; i < N_; i++) {
          V[j] += wts[i] * sigma[offset+i];
        }
        offset += N_;
      }
    }
  }

  template <class Real, Integer Order> const std::string& DiscElastance<Real,Order>::Name() const {
    const double R = (double)this->disc_panels.DiscRadius();
    static auto name = std::string("LaplaceElastance") + std::to_string(Order) + "-R" + std::to_string(R);
    return name;
  }

  template <class Real, Integer Order> void DiscElastance<Real,Order>::BuildInteracBlock(Matrix<Real>& M, const DiscPanelLst<Real,Order> panel_lst, const typename DiscPanelLst<Real,Order>::NearData& interac_block, const Real tol) const {
    const Long disc_idx0 = interac_block.disc_idx0;
    const Long disc_idx1 = interac_block.disc_idx1;
    const Long disc_panel_start0 = interac_block.panel_idx_range0[0];
    const Long disc_panel_start1 = interac_block.panel_idx_range1[0];
    const Long disc_panel_cnt0 = interac_block.panel_idx_range0[1] - disc_panel_start0;
    const Long disc_panel_cnt1 = interac_block.panel_idx_range1[1] - disc_panel_start1;
    const Long N0 = disc_panel_cnt0 * Order;
    const Long N1 = disc_panel_cnt1 * Order;

    Vector<Real> Xt((N0 + N1) * COORD_DIM);
    { // Set Xt
      const auto& X0 = panel_lst.SurfCoord(disc_idx0);
      const auto& X1 = panel_lst.SurfCoord(disc_idx1);
      for (Long i = 0; i < N0*COORD_DIM; i++) Xt[             i] = X0[disc_panel_start0 * Order*COORD_DIM+i];
      for (Long i = 0; i < N1*COORD_DIM; i++) Xt[N0*COORD_DIM+i] = X1[disc_panel_start1 * Order*COORD_DIM+i];
    }

    Matrix<Real> Md(N0+N1, N0+N1);
    { // Set Md
      static_assert(Laplace2D_DxU::SrcDim() == 1, "");
      static_assert(Laplace2D_DxU::TrgDim() == 1, "");
      Matrix<Real> Md0(N0, N0+N1, Md[0], false);
      Matrix<Real> Md1(N1, N0+N1, Md[N0], false);
      const Long panel_start0 = panel_lst.PanelIdxOffset(disc_idx0) + disc_panel_start0;
      const Long panel_start1 = panel_lst.PanelIdxOffset(disc_idx1) + disc_panel_start1;
      panel_lst.template LayerPotentialMatrix<Laplace2D_DxU,-1>(Md0, Xt, tol, panel_start0, panel_start0 + disc_panel_cnt0);
      panel_lst.template LayerPotentialMatrix<Laplace2D_DxU,-1>(Md1, Xt, tol, panel_start1, panel_start1 + disc_panel_cnt1);
      for (Long i = 0; i < Md.Dim(0); i++) Md[i][i] += 0.5;
    }

    Matrix<Real> VVt(N0+N1, N0+N1);
    { // Set VVt
      const auto& W0 = panel_lst.SurfWts(disc_idx0);
      const auto& W1 = panel_lst.SurfWts(disc_idx1);
      const Long offset0 = disc_panel_start0 * Order;
      const Long offset1 = disc_panel_start1 * Order;

      VVt.SetZero();
      for (Long i = 0; i < N0; i++) {
        for (Long j = 0; j < N0; j++) {
          VVt[i][j] += W0[offset0+i];
        }
      }
      for (Long i = 0; i < N1; i++) {
        for (Long j = 0; j < N1; j++) {
          VVt[N0+i][N0+j] = W1[offset1+i];
        }
      }
    }

    M = Md + VVt;
  }

  template <class Real, Integer Order> void DiscElastance<Real,Order>::ApplyBIOpDirect(Vector<Real>* U, const Vector<Real>& sigma) const {
    Vector<Real> sigma_near, sigma_far;
    Vector<Real> sigma_near_, sigma_far_;
    this->Split(&sigma_near, &sigma_far, sigma);
    this->Merge(&sigma_near_, sigma_near, Vector<Real>());
    this->Merge(&sigma_far_, Vector<Real>(), sigma_far);

    Vector<Real> Udl;
    LaplaceDL_BIOp_far .ComputePotential(Udl, sigma_far );
    if (sigma_near.Dim()) {
      Vector<Real> Udl_near, Udl_near_;
      LaplaceDL_BIOp_near.ComputePotential(Udl_near, sigma_near);
      this->Merge(&Udl_near_, Vector<Real>(), Udl_near);
      Udl += Udl_near_;
    }
    (*U) = 0.5*sigma_far_ + Udl;

    Long offset = 0;
    Vector<Real> VVt_near(U->Dim());
    Vector<Real> VVt_far (U->Dim());
    for (Long k = 0; k < this->disc_panels.DiscCount(); k++) {
      const auto& wts = this->disc_panels.SurfWts(k);
      const Long N = wts.Dim();

      Real Qn = 0, Qf = 0;
      if (sigma_near_.Dim()) {
        for (Long i = 0; i < N; i++) Qn += wts[i] * sigma_near_[offset+i];
      }
      for (Long i = 0; i < N; i++) Qf += wts[i] * sigma_far_ [offset+i];

      for (Long i = 0; i < N; i++) VVt_near[offset+i] = Qn;
      for (Long i = 0; i < N; i++) VVt_far [offset+i] = Qf;
      offset += N;
    }
    { // Set near-near block of VVt_near to zero
      Vector<Real> Vtmp;
      this->Split(nullptr, &Vtmp, VVt_near);
      this->Merge(&VVt_near, Vector<Real>(), Vtmp);
    }
    (*U) += VVt_near + VVt_far;
  }

}
