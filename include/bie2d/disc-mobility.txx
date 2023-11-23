namespace sctl {

  template <class Real, Integer Order> DiscMobility<Real,Order>::DiscMobility(const Comm& comm_) : ICIP_Base(comm_), StokesSL_BIOp(this->comm), StokesDL_BIOp(this->comm) {}

  template <class Real, Integer Order> void DiscMobility<Real,Order>::Init(const Vector<Real>& Xc, const Real R, const Real tol, const ICIPType icip_type) {
    ICIP_Base::Init(Xc, R, tol, icip_type);
    RigidVelocityBasis(V0, this->disc_panels);

    const bool exclude_near = (icip_type == ICIPType::Compress || icip_type == ICIPType::Precond);
    StokesDL_BIOp.Init(this->disc_panels, exclude_near, tol);
    StokesSL_BIOp.Init(this->disc_panels, false, tol);
  }

  template <class Real, Integer Order> void DiscMobility<Real,Order>::Solve(Vector<Real>& V, const Vector<Real>& F, const Vector<Real>& Vs, const Real gmres_tol, const Long gmres_max_iter) {
    const Long Ndisc = this->disc_panels.DiscCount();
    const Long N = this->disc_panels.Size() * Order * COORD_DIM;
    const Real R = this->disc_panels.DiscRadius();

    Vector<Real> nu; // boundary force
    { // Set nu
      nu.ReInit(N);
      Long offset = 0;
      for (Long i = 0; i < Ndisc; i++) {
        const Long N_ = this->disc_panels.SurfWts(i).Dim()*COORD_DIM;
        const Real inv_circumf = 1/(2*const_pi<Real>()*R);
        for (Long j = 0; j < N_; j++) {
          nu[offset+j] = 0;
          for (Integer k = 0; k < 3; k++) {
            nu[offset+j] += F[i*3+k] * V0[k][offset+j] * inv_circumf;
          }
        }
        offset += N_;
      }
    }

    Vector<Real> U0; // completion flow velocity
    StokesSL_BIOp.ComputePotential(U0, nu);

    Vector<Real> sigma; // unknown density
    this->SolveBIE(sigma, U0, gmres_tol, gmres_max_iter);

    if (V.Dim() != Ndisc*3) V.ReInit(Ndisc*3);
    { // recover translation and rotation from sigma
      V = 0;
      Long offset = 0;
      for (Long j = 0; j < Ndisc; j++) {
        const auto& wts = this->disc_panels.SurfWts(j);
        const Long N_ = wts.Dim() * COORD_DIM;
        const Real inv_circumf = 1/(2*const_pi<Real>()*R);
        for (Long i = 0; i < N_; i++) {
          const Long ii = (i>>1);
          for (Integer k = 0; k < 3; k++) {
            V[j*3+k] += inv_circumf * wts[ii] * V0[k][offset+i] * sigma[offset+i];
          }
        }
        offset += N_;
      }
    }
  }

  template <class Real, Integer Order> const std::string& DiscMobility<Real,Order>::Name() const {
    const double R = (double)this->disc_panels.DiscRadius();
    static auto name = std::string("StokesMobility") + std::to_string(Order) + "-R" + std::to_string(R);
    return name;
  }

  template <class Real, Integer Order> void DiscMobility<Real,Order>::BuildInteracBlock(Matrix<Real>& M, const DiscPanelLst<Real,Order>& panel_lst, const typename DiscPanelLst<Real,Order>::NearData& interac_block, const Real tol) const {
    const Long disc_idx0 = interac_block.disc_idx0;
    const Long disc_idx1 = interac_block.disc_idx1;
    const Long disc_panel_start0 = interac_block.panel_idx_range0[0];
    const Long disc_panel_start1 = interac_block.panel_idx_range1[0];
    const Long disc_panel_cnt0 = interac_block.panel_idx_range0[1] - disc_panel_start0;
    const Long disc_panel_cnt1 = interac_block.panel_idx_range1[1] - disc_panel_start1;
    const Long N0 = disc_panel_cnt0 * Order;
    const Long N1 = disc_panel_cnt1 * Order;

    Matrix<Real> V0, V1;
    RigidVelocityBasis(V0, panel_lst, disc_idx0);
    RigidVelocityBasis(V1, panel_lst, disc_idx1);

    Vector<Real> Xt((N0 + N1) * COORD_DIM);
    { // Set Xt
      const auto& X0 = panel_lst.SurfCoord(disc_idx0);
      const auto& X1 = panel_lst.SurfCoord(disc_idx1);
      for (Long i = 0; i < N0*COORD_DIM; i++) Xt[             i] = X0[disc_panel_start0 * Order*COORD_DIM+i];
      for (Long i = 0; i < N1*COORD_DIM; i++) Xt[N0*COORD_DIM+i] = X1[disc_panel_start1 * Order*COORD_DIM+i];
    }

    static constexpr Integer DOF = 2;
    Matrix<Real> Md((N0+N1)*DOF, (N0+N1)*DOF);
    { // Set Md
      static_assert(Stokes2D_DxU::SrcDim() == DOF, "");
      static_assert(Stokes2D_DxU::TrgDim() == DOF, "");
      Matrix<Real> Md0(N0*DOF, (N0+N1)*DOF, Md[0], false);
      Matrix<Real> Md1(N1*DOF, (N0+N1)*DOF, Md[N0*DOF], false);
      const Long panel_start0 = panel_lst.PanelIdxOffset(disc_idx0) + disc_panel_start0;
      const Long panel_start1 = panel_lst.PanelIdxOffset(disc_idx1) + disc_panel_start1;
      panel_lst.template LayerPotentialMatrix<Stokes2D_DxU,-1>(Md0, Xt, tol, panel_start0, panel_start0 + disc_panel_cnt0);
      panel_lst.template LayerPotentialMatrix<Stokes2D_DxU,-1>(Md1, Xt, tol, panel_start1, panel_start1 + disc_panel_cnt1);
      for (Long i = 0; i < Md.Dim(0); i++) Md[i][i] += 0.5;
    }

    Matrix<Real> VVt((N0+N1)*DOF, (N0+N1)*DOF);
    { // Set VVt
      const Real inv_circumf = 1/(2*const_pi<Real>()*panel_lst.DiscRadius());
      const auto& W0 = panel_lst.SurfWts(disc_idx0);
      const auto& W1 = panel_lst.SurfWts(disc_idx1);
      const Long offset0 = disc_panel_start0 * Order;
      const Long offset1 = disc_panel_start1 * Order;

      VVt.SetZero();
      for (Long i = 0; i < N0*DOF; i++) {
        const Real wt = W0[offset0+i/DOF] * inv_circumf;
        for (Long j = 0; j < N0*DOF; j++) {
          for (Long k = 0; k < 3; k++) {
            VVt[i][j] += V0[k][offset0*DOF+i] * V0[k][offset0*DOF+j] * wt;
          }
        }
      }
      for (Long i = 0; i < N1*DOF; i++) {
        const Real wt = W1[offset1+i/DOF] * inv_circumf;
        for (Long j = 0; j < N1*DOF; j++) {
          for (Long k = 0; k < 3; k++) {
            VVt[N0*DOF+i][N0*DOF+j] += V1[k][offset1*DOF+i] * V1[k][offset1*DOF+j] * wt;
          }
        }
      }
    }

    M = Md + VVt;
  }

  template <class Real, Integer Order> void DiscMobility<Real,Order>::ApplyBIOpDirect(Vector<Real>* U, const Vector<Real>& sigma) const {
    Vector<Real> sigma_near, sigma_far;
    Vector<Real> sigma_near_, sigma_far_;
    this->Split(&sigma_near, &sigma_far, sigma);
    this->Merge(&sigma_near_, sigma_near, Vector<Real>());
    this->Merge(&sigma_far_, Vector<Real>(), sigma_far);

    Vector<Real> Udl;
    StokesDL_BIOp.ComputePotential(Udl, sigma);
    (*U) = 0.5*sigma_far_ + Udl;

    const auto apply_wts = [](Vector<Real>& sigma, const Vector<Real>& wt, const Real scal) {
      const Long N = wt.Dim();
      const Long dof = sigma.Dim() / N;
      SCTL_ASSERT(sigma.Dim() == N * dof);
      for (Long i = 0; i < N; i++) {
        for (Long k = 0; k < dof; k++) {
          sigma[i*dof+k] *= wt[i] * scal;
        }
      }
    };
    const Real inv_circumf = 1/(2*const_pi<Real>()*this->disc_panels.DiscRadius());
    apply_wts(sigma_far_, this->disc_panels.SurfWts(-1), inv_circumf);
    apply_wts(sigma_near_, this->disc_panels.SurfWts(-1), inv_circumf);
    this->disc_wise_outer_product(*U, sigma_far_, sigma_near_, V0, V0);
  }

  template <class Real, Integer Order> void DiscMobility<Real,Order>::RigidVelocityBasis(Matrix<Real>& V, const DiscPanelLst<Real,Order>& disc_panels, const Long disc_idx) {
    const Long N = (disc_idx < 0 ? disc_panels.Size()*Order : disc_panels.SurfWts(disc_idx).Dim());
    if (V.Dim(0) != 3 || V.Dim(1) != N*COORD_DIM) V.ReInit(3, N*COORD_DIM);

    Long offset = 0;
    const Long disc_idx0 = (disc_idx < 0 ? 0 : disc_idx);
    const Long Ndisc = (disc_idx < 0 ? disc_panels.DiscCount() : 1);
    for (Long k = disc_idx0; k < disc_idx0 + Ndisc; k++) {
      const Real xx = std::get<0>(disc_panels.DiscCoord(k));
      const Real yy = std::get<1>(disc_panels.DiscCoord(k));

      const auto& X = disc_panels.SurfCoord(k);
      const Long N_ = X.Dim() / COORD_DIM;

      const Real Rinv = 1/disc_panels.DiscRadius();
      for (Long i = 0; i < N_; i++) {
        // Translation along X
        V[0][(offset+i)*COORD_DIM+0] = 1;
        V[0][(offset+i)*COORD_DIM+1] = 0;

        // Translation along Y
        V[1][(offset+i)*COORD_DIM+0] = 0;
        V[1][(offset+i)*COORD_DIM+1] = 1;

        // Rotation
        V[2][(offset+i)*COORD_DIM+0] =-(X[i*COORD_DIM+1]-yy) * Rinv;
        V[2][(offset+i)*COORD_DIM+1] = (X[i*COORD_DIM+0]-xx) * Rinv;
      }
      offset += N_;
    }
  }

}
