namespace sctl {

  template <class Real, Integer Order> DiscCapacitance<Real,Order>::DiscCapacitance(const Comm& comm_) : ICIP_Base(comm_), LaplaceSL_BIOp(this->comm), LaplaceDL_BIOp(this->comm) {}

  template <class Real, Integer Order> void DiscCapacitance<Real,Order>::Init(const Vector<Real>& Xc, const Real R, const Real tol, const ICIPType icip_type) {
    ICIP_Base::Init(Xc, R, tol, icip_type);

    const bool exclude_near = (icip_type == ICIPType::Compress || icip_type == ICIPType::Precond);
    LaplaceSL_BIOp.Init(this->disc_panels, exclude_near, tol);
    LaplaceDL_BIOp.Init(this->disc_panels, exclude_near, tol);
  }

  template <class Real, Integer Order> void DiscCapacitance<Real,Order>::Solve(Vector<Real>& Q, const Vector<Real>& V, const Real gmres_tol, const Long gmres_max_iter) {
    const Long Ndisc = this->disc_panels.DiscCount();
    const Long N = this->disc_panels.Size() * Order;

    Vector<Real> v; // pointwise potential
    { // Set v
      v.ReInit(N);
      Long offset = 0;
      for (Long i = 0; i < Ndisc; i++) {
        const Long N_ = this->disc_panels.SurfWts(i).Dim();
        for (Long j = 0; j < N_; j++) v[offset+j] = V[i];
        offset += N_;
      }
    }

    Vector<Real> sigma; // unknown density
    this->SolveBIE(sigma, v, gmres_tol, gmres_max_iter);

    if (Q.Dim() != Ndisc) Q.ReInit(Ndisc);
    { // get charge on each disc
      Q = 0;
      Long offset = 0;
      for (Long j = 0; j < Ndisc; j++) {
        const auto& wts = this->disc_panels.SurfWts(j);
        const Long N_ = wts.Dim();
        for (Long i = 0; i < N_; i++) {
          Q[j] += wts[i] * sigma[offset+i];
        }
        offset += N_;
      }
    }
  }

  template <class Real, Integer Order> const std::string& DiscCapacitance<Real,Order>::Name() const {
    const double R = (double)this->disc_panels.DiscRadius();
    static auto name = std::string("LaplaceCapacitance") + std::to_string(Order) + "-R" + std::to_string((double)R);
    return name;
  }

  template <class Real, Integer Order> void DiscCapacitance<Real,Order>::BuildInteracBlock(Matrix<Real>& M, const DiscPanelLst<Real,Order>& panel_lst, const typename DiscPanelLst<Real,Order>::NearData& interac_block, const Real tol) const {
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
    Matrix<Real> Ms(N0+N1, N0+N1);
    { // Set Md, Ms
      static_assert(Laplace2D_DxU::SrcDim() == 1, "");
      static_assert(Laplace2D_DxU::TrgDim() == 1, "");
      const Long panel_start0 = panel_lst.PanelIdxOffset(disc_idx0) + disc_panel_start0;
      const Long panel_start1 = panel_lst.PanelIdxOffset(disc_idx1) + disc_panel_start1;

      Matrix<Real> Md0(N0, N0+N1, Md[0], false);
      Matrix<Real> Md1(N1, N0+N1, Md[N0], false);
      panel_lst.template LayerPotentialMatrix<Laplace2D_DxU,-1>(Md0, Xt, tol, panel_start0, panel_start0 + disc_panel_cnt0);
      panel_lst.template LayerPotentialMatrix<Laplace2D_DxU,-1>(Md1, Xt, tol, panel_start1, panel_start1 + disc_panel_cnt1);
      for (Long i = 0; i < Md.Dim(0); i++) Md[i][i] += 0.5;

      Matrix<Real> Ms0(N0, N0+N1, Ms[0], false);
      Matrix<Real> Ms1(N1, N0+N1, Ms[N0], false);
      panel_lst.template LayerPotentialMatrix<Laplace2D_FxU,-1>(Ms0, Xt, tol, panel_start0, panel_start0 + disc_panel_cnt0);
      panel_lst.template LayerPotentialMatrix<Laplace2D_FxU,-1>(Ms1, Xt, tol, panel_start1, panel_start1 + disc_panel_cnt1);
    }

    M = Md + Ms;
  }

  template <class Real, Integer Order> void DiscCapacitance<Real,Order>::ApplyBIOpDirect(Vector<Real>* U, const Vector<Real>& sigma) const {
    Vector<Real> sigma_near, sigma_far, sigma_far_;
    this->Split(&sigma_near, &sigma_far, sigma);
    this->Merge(&sigma_far_, Vector<Real>(), sigma_far);

    Vector<Real> Udl, Usl;
    LaplaceDL_BIOp.ComputePotential(Udl, sigma);
    LaplaceSL_BIOp.ComputePotential(Usl, sigma);
    (*U) = 0.5*sigma_far_ + Udl + Usl;
  }

}
