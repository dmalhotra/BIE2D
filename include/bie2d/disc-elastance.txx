namespace sctl {

  template <class Real, Integer Order> DiscElastance<Real,Order>::DiscElastance(const Comm& comm_) : ICIP_Base(comm_), LaplaceSL_BIOp(this->comm), LaplaceDL_BIOp(this->comm) {}

  template <class Real, Integer Order> void DiscElastance<Real,Order>::Init(const Vector<Real>& Xc, const Real R, const Real tol, const ICIPType icip_type) {
    ICIP_Base::Init(Xc, R, tol, icip_type);

    const bool exclude_near = (icip_type == ICIPType::Compress || icip_type == ICIPType::Precond);
    LaplaceDL_BIOp.Init(this->disc_panels, exclude_near, tol);
    LaplaceSL_BIOp.Init(this->disc_panels, false, tol);
  }

  template <class Real, Integer Order> void DiscElastance<Real,Order>::Solve(Vector<Real>& V, const Vector<Real>& Q, const Real gmres_tol, const Long gmres_max_iter) {
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

    Vector<Real> sigma; // unknown density
    this->SolveBIE(sigma, U0, gmres_tol, gmres_max_iter);

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

  template <class Real, Integer Order> void DiscElastance<Real,Order>::BuildInteracBlock(Matrix<Real>& M, const DiscPanelLst<Real,Order>& panel_lst, const typename DiscPanelLst<Real,Order>::NearData& interac_block, const Real tol) const {
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
    LaplaceDL_BIOp.ComputePotential(Udl, sigma);
    (*U) = 0.5*sigma_far_ + Udl;

    Matrix<Real> one(1,sigma.Dim()); one = 1;
    this->disc_wise_outer_product(*U, sigma_far_*this->disc_panels.SurfWts(-1), (sigma_near_.Dim() ? sigma_near_*this->disc_panels.SurfWts(-1) : Vector<Real>()), one, one);
  }

}
