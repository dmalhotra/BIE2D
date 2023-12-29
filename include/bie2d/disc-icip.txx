namespace sctl {

  template <class Real, Integer Order> ICIP<Real,Order>::ICIP(const Comm& comm_) : comm(comm_), solver(comm, true) {}

  template <class Real, Integer Order> ICIP<Real,Order>::~ICIP() {}

  template <class Real, Integer Order> void ICIP<Real,Order>::Init(const Vector<Real>& Xc, const Real R, const Real tol, const ICIPType icip_type) {
    disc_panels.Init(Xc, R, (icip_type==ICIPType::Adaptive), d_max);
    icip_type_ = icip_type;
    Kcorrec.ReInit(0);
    Rprecon.ReInit(0);
    tol_ = tol;

    { // setup split into near and far panels
      X = disc_panels.SurfCoord(-1);
      Split(&Xnear, &Xfar, X);
      panels_near.Init(Xnear);
      panels_far.Init(Xfar);
    }

    { // Set sqrt_wts, rsqrt_wts
      const auto& wts = disc_panels.SurfWts(-1);
      const Long N = wts.Dim();
      if ( sqrt_wts.Dim() != N)  sqrt_wts.ReInit(N);
      if (rsqrt_wts.Dim() != N) rsqrt_wts.ReInit(N);
      for (Long i = 0; i < N; i++) {
        const Real sqrt_w = sqrt<Real>(wts[i]);
        sqrt_wts[i] = sqrt_w;
        rsqrt_wts[i] = 1/sqrt_w;
      }
    }
  }

  template <class Real, Integer Order> const DiscPanelLst<Real,Order>& ICIP<Real,Order>::GetPanelList() const {
    return disc_panels;
  }

  template <class Real, Integer Order> void ICIP<Real,Order>::BuildCompression(Matrix<Real>* R, Matrix<Real>* Rinv, const Real x0, const Real y0, const Real x1, const Real y1, const Real radius, const Real tol) const {

    DiscPanelLst<Real,Order> panels_coarse, panels_fine;
    const auto Xc = Vector<Real>({x0,y0,x1,y1});
    panels_coarse.Init(Xc, radius, false, d_max);
    panels_fine  .Init(Xc, radius,  true, d_max);

    if (panels_coarse.GetNearList().Dim() && panels_fine.GetNearList().Dim()) {
      SCTL_ASSERT(panels_coarse.GetNearList().Dim() == 1);
      SCTL_ASSERT(panels_fine  .GetNearList().Dim() == 1);
      auto interac_coarse = panels_coarse.GetNearList()[0];
      auto interac_fine   = panels_fine  .GetNearList()[0];
      if (interac_coarse.disc_idx0 > interac_coarse.disc_idx1) {
        std::swap(interac_coarse.disc_idx0, interac_coarse.disc_idx1);
        std::swap(interac_coarse.panel_idx_range0[0], interac_coarse.panel_idx_range1[0]);
        std::swap(interac_coarse.panel_idx_range0[1], interac_coarse.panel_idx_range1[1]);
      }
      if (interac_fine  .disc_idx0 > interac_fine  .disc_idx1) {
        std::swap(interac_fine.disc_idx0, interac_fine.disc_idx1);
        std::swap(interac_fine.panel_idx_range0[0], interac_fine.panel_idx_range1[0]);
        std::swap(interac_fine.panel_idx_range0[1], interac_fine.panel_idx_range1[1]);
      }

      Matrix<Real> Kf, Kf_inv;
      BuildInteracBlock(Kf, panels_fine, interac_fine, tol);
      Kf_inv = Matrix<Real>(Kf).pinv(tol);
      SCTL_ASSERT(Kf.Dim(0) == Kf.Dim(1));

      Matrix<Real> Wc_inv, Wf;
      { // Build Wc_inv, Wf
        const auto get_wts = [](Vector<Real>& W, const DiscPanelLst<Real,Order>& panels, const typename DiscPanelLst<Real,Order>::NearData& interac) {
          const Long offset0 = interac.panel_idx_range0[0]*Order;
          const Long offset1 = interac.panel_idx_range1[0]*Order;
          const Long count0 = interac.panel_idx_range0[1]*Order - offset0;
          const Long count1 = interac.panel_idx_range1[1]*Order - offset1;
          if (W.Dim() != count0+count1) W.ReInit(count0+count1);

          const auto& W0 = panels.SurfWts(interac.disc_idx0);
          const auto& W1 = panels.SurfWts(interac.disc_idx1);
          for (Long i = 0; i < count0; i++) W[       i] = W0[offset0+i];
          for (Long i = 0; i < count1; i++) W[count0+i] = W1[offset1+i];
        };
        Vector<Real> Wc_, Wf_;
        get_wts(Wc_, panels_coarse, interac_coarse);
        get_wts(Wf_, panels_fine  , interac_fine  );
        const Long Nc = Wc_.Dim();
        const Long Nf = Wf_.Dim();

        const Integer dof = Kf.Dim(0) / Nf;
        SCTL_ASSERT(Kf.Dim(0) == Nf * dof);

        Wc_inv.ReInit(Nc*dof, Nc*dof); Wc_inv.SetZero();
        for (Long i = 0; i < Nc; i++) {
          for (Long k = 0; k < dof; k++) {
            Wc_inv[i*dof+k][i*dof+k] = 1/Wc_[i];
          }
        }

        Wf.ReInit(Nf*dof, Nf*dof); Wf.SetZero();
        for (Long i = 0; i < Nf; i++) {
          for (Long k = 0; k < dof; k++) {
            Wf[i*dof+k][i*dof+k] = Wf_[i];
          }
        }
      }

      Matrix<Real> P;
      { // Build P
        const auto get_theta_breaks = [](Vector<Real>& T0, Vector<Real>& T1, const DiscPanelLst<Real,Order>& panels, const typename DiscPanelLst<Real,Order>::NearData& interac) {
          const Long offset0 = interac.panel_idx_range0[0];
          const Long offset1 = interac.panel_idx_range1[0];
          const Long count0 = interac.panel_idx_range0[1] - offset0;
          const Long count1 = interac.panel_idx_range1[1] - offset1;
          if (T0.Dim() != count0+1) T0.ReInit(count0+1);
          if (T1.Dim() != count1+1) T1.ReInit(count1+1);

          const auto& TT0 = panels.SurfTheta(interac.disc_idx0);
          const auto& TT1 = panels.SurfTheta(interac.disc_idx1);
          for (Long i = 0; i < count0; i++) T0[i] = TT0[offset0+i];
          for (Long i = 0; i < count1; i++) T1[i] = TT1[offset1+i];
          T0[count0] = (offset0+count0<TT0.Dim() ? TT0[offset0+count0] : TT0[0]+2*const_pi<Real>());
          T1[count1] = (offset1+count1<TT1.Dim() ? TT1[offset1+count1] : TT1[0]+2*const_pi<Real>());
        };
        Vector<Real> Tc0, Tc1, Tf0, Tf1;
        get_theta_breaks(Tc0, Tc1, panels_coarse, interac_coarse);
        get_theta_breaks(Tf0, Tf1, panels_fine  , interac_fine  );

        const auto& nds = panels_coarse.PanelNds();
        const auto build_interp_mat = [&nds](Matrix<Real>& P, Vector<Real> Tc, Vector<Real> Tf) {
          const Long Nc = (Tc.Dim()-1);
          const Long Nf = (Tf.Dim()-1);
          SCTL_ASSERT(Nc && Nf);

          while (Tc[0] - Tf[0] > const_pi<Real>()) Tc -= const_pi<Real>();
          while (Tf[0] - Tc[0] > const_pi<Real>()) Tf -= const_pi<Real>();

          if (P.Dim(0) != Nc*Order && P.Dim(1) != Nf*Order) P.ReInit(Nc*Order, Nf*Order);
          P.SetZero();

          Long src_panel_idx = 0;
          Vector<Real> nds_src, nds_trg, P_(Order*Order); // temporary storage
          for (Long trg_panel_idx = 0; trg_panel_idx < Nf; trg_panel_idx++) {
            while (src_panel_idx<Nc-1 && Tc[src_panel_idx+1]<Tf[trg_panel_idx+1]) src_panel_idx++;
            SCTL_ASSERT(Tc[src_panel_idx] < Tf[trg_panel_idx+1]);
            SCTL_ASSERT(Tf[trg_panel_idx] < Tc[src_panel_idx+1]);

            nds_trg = nds;
            nds_trg *= (Tf[trg_panel_idx+1]-Tf[trg_panel_idx]);
            nds_trg += Tf[trg_panel_idx];

            nds_src = nds;
            nds_src *= (Tc[src_panel_idx+1]-Tc[src_panel_idx]);
            nds_src += Tc[src_panel_idx];

            LagrangeInterp<Real>::Interpolate(P_, nds_src, nds_trg);
            for (Long i = 0; i < Order; i++) {
              for (Long j = 0; j < Order; j++) {
                P[src_panel_idx*Order+i][trg_panel_idx*Order+j] = P_[i*Order+j];
              }
            }
          }
        };
        Matrix<Real> P0, P1;
        build_interp_mat(P0, Tc0, Tf0);
        build_interp_mat(P1, Tc1, Tf1);

        const Long Nc0 = P0.Dim(0);
        const Long Nc1 = P1.Dim(0);
        const Long Nf0 = P0.Dim(1);
        const Long Nf1 = P1.Dim(1);
        const Long dof = Kf.Dim(0) / (Nf0+Nf1);
        SCTL_ASSERT(Kf.Dim(0) == (Nf0+Nf1)*dof);

        P.ReInit((Nc0+Nc1)*dof, (Nf0+Nf1)*dof); P.SetZero();
        for (Long i = 0; i < Nc0; i++) {
          for (Long j = 0; j < Nf0; j++) {
            for (Long k = 0; k < dof; k++) {
              P[i*dof+k][j*dof+k] = P0[i][j];
            }
          }
        }
        for (Long i = 0; i < Nc1; i++) {
          for (Long j = 0; j < Nf1; j++) {
            for (Long k = 0; k < dof; k++) {
              P[(Nc0+i)*dof+k][(Nf0+j)*dof+k] = P1[i][j];
            }
          }
        }
      }
      
      const auto R_ = P * Kf_inv * Wf * P.Transpose() * Wc_inv;

      if (R != nullptr) (*R) = R_;
      if (Rinv != nullptr) (*Rinv) = Matrix<Real>(R_).pinv(tol);
    }
  }

  template <class Real, Integer Order> void ICIP<Real,Order>::GetPrecondBlock(Matrix<Real>* R, Matrix<Real>* Rinv, const Real x0, const Real y0, const Real x1, const Real y1, const Real radius) const {
    Real tol = machine_eps<Real>();
    char sep = '_'; // How files are separated according to the Legendre node number

    // @Mariana TODO: get R and Rinv from interpolation instead of computing on-the-fly.
    // Don't know if this works - still testing (problems with C++)
    //
    // IF THERE IS NO FILE THEN BUILD THE PRECOMPUTED R and Rinv
    const auto& leg_nds = LegQuadRule<Real>::ComputeNds(InterpOrder); // compute Legendre nodes in [-1,1]
    const char *fRnum, *fRinvnum;
    
    std::string fR = "include/bie2d/precomputed/R";
    std::string fRinv = "include/bie2d/precomputed/Rinv";
    std::string s1, s2, s3, s4;
    if( icip_type_ == ICIPType::Precond){
      s1 = fR + sep + std::to_string(0) + ".bin";
      s2 = fR + sep + std::to_string(InterpOrder-1) + ".bin";
    }
    else{
      s3 = fRinv + sep + std::to_string(0) + ".bin";
      s4 = fRinv + sep + std::to_string(InterpOrder-1) + ".bin";
    }
    // CHECK IF NECESSARY FILES EXIST
    if( (access( s1.c_str(), F_OK ) == -1 || access( s2.c_str(), F_OK ) == -1) && R != nullptr ){
      // IF THERE IS NO FIRST FILE OR LAST FILE THEN BUILD R (COMPRESS)
      std::cout << "\nCreating files for InterpOrder: "<< InterpOrder << "\n";
      Matrix<Real> R0[InterpOrder], Rinv0[InterpOrder];
      for (Long i = 0; i < InterpOrder; i++) { // loop over interpolation nodes
       // get name of file
       s2 = fR + sep + std::to_string(i) + ".bin";
       fRnum = s2.c_str();
       // transform to log Legendre and get the new centers of the discs
       const Real x = radius + radius*(d_min/2) * exp(log(d_max/d_min)*(leg_nds[i] + 1)/2);
       //
       // Compute R and for this distance in particular
       std::cout << "\nBuilding precompression for log Legendre for i: "<< i << "\n";
       std::cout << "Size of R0: " << sizeof(R0[i]) << " size of Rinv0: " << sizeof(Rinv0[i]) << "\n";
       std::cout << "x: " << x << " radius: " << radius << " tol: " << tol;
       BuildCompression(&R0[i], &Rinv0[i], -x, 0, x, 0, radius, tol);
     
       // Save as a binary file using the SCTL library
       std::cout << "Save to file: \n";
       std::cout << fRnum << std::endl;
       R0[i].template Write<double>(fRnum); // to specify saving in double-precision format

      }
    }
    if( (access( s3.c_str(), F_OK ) == -1 || access( s4.c_str(), F_OK ) == -1) && Rinv != nullptr ){
      // IF THERE IS NO FIRST FILE OR LAST FILE THEN BUILD Rinv (Preconditioning)
      Matrix<Real> R0[InterpOrder], Rinv0[InterpOrder];
      for (Long i = 0; i < InterpOrder; i++) { // loop over interpolation nodes
       // get name of file
       s4 = fRinv + sep + std::to_string(i) + ".bin";
       fRinvnum = s4.c_str();
       // transform to log Legendre and get the new centers of the discs
       const Real x = radius + radius*(d_min/2) * exp(log(d_max/d_min)*(leg_nds[i] + 1)/2);
       //
       // Compute R and Rinv for this distance in particular
       std::cout << "\nBuilding precompression for log Legendre\n";
       BuildCompression(&R0[i], &Rinv0[i], -x, 0, x, 0, radius, tol);
     
       // Save as a binary file using the SCTL library
       std::cout << "Save to file: \n";
       std::cout << fRinvnum << std::endl;
       Rinv0[i].template Write<double>(fRinvnum);
      }
    }
    
    // INTERPOLATE
    const Real d = log((sqrt( (x0-x1)*(x0-x1 ) + (y0-y1)*(y0-y1) ) - 2*radius)/(radius*d_min))/log(d_max/d_min);
    
    // Get w for Barycentric interpolation
    Vector<Real> w(InterpOrder);
    Real maxNorm = 0;
    for( Integer j = 0; j<InterpOrder; j++){
      Real w_inv = 1;
      Real this_node = leg_nds[j];
      for( Integer k = 0  ; k< j         ; k++) w_inv *= leg_nds[k] - this_node;
      for( Integer k = j+1; k<InterpOrder; k++) w_inv *= leg_nds[k] - this_node;
      w[j] = 1/w_inv;
      if( w[j] > maxNorm) maxNorm = w[j];
    }
    w = w/maxNorm;
    //std::cout << "maxNorm: " << maxNorm << "\n";
    
    // Build R interpolated if it's not a null pointer
    Matrix<Real> Rtemp;
    Matrix<Real> RinvTemp;
    if( R!= nullptr ){
      Matrix<Real> Ri[InterpOrder];
      std::cout << "Interpolating R for " << icip_type_ << "\n";
      s1 = fR + sep + std::to_string(0) + ".bin";
      fRnum = s1.c_str();
      Ri[0].template Read<double>(fRnum); // First matrix
      Long init = 1;
      if( fabs(d - leg_nds[0])<tol){
        init = InterpOrder;
	Rtemp = Ri[0];
      }
      Real coef = w[0]/(d - leg_nds[0]);
      Real denom = w[0]/(d - leg_nds[0]);
      Rtemp = coef*Ri[0];
      // Iterate
      for( Long i = init; i<InterpOrder; i++){
        s1 = fR + sep + std::to_string(i) + ".bin";
	fRnum = s1.c_str();
        Ri[i].template Read<double>(fRnum); // Read matrix from file
        if( fabs(d - leg_nds[i])<tol){
	  Rtemp = Ri[i];
	  break;
	}
        coef = w[i]/(d - leg_nds[i]);
	denom += coef;
	Rtemp += coef*Ri[i];
      }
      Rtemp = Rtemp/denom;
    }
    // Build Rinv interpolated if it's not a null pointer
    if( Rinv!= nullptr ){
      Matrix<Real> Rinvi[InterpOrder];
      std::cout << "Interpolating Rinv for " << icip_type_ << "\n";
      s3 = fRinv + sep + std::to_string(0) + ".bin";
      fRnum = s3.c_str();
      Rinvi[0].template Read<double>(fRnum); // First matrix
      Long init = 1;
      if( fabs(d - leg_nds[0])<tol){
        init = InterpOrder;
	RinvTemp = Rinvi[0];
      }
      Real coef = w[0]/(d - leg_nds[0]);
      Real denom = w[0]/(d - leg_nds[0]);
      RinvTemp = coef*Rinvi[0];
      // Iterate
      for( Long i = init; i<InterpOrder; i++){
        s3 = fRinv + sep + std::to_string(i) + ".bin";
	fRinvnum = s3.c_str();
        Rinvi[i].template Read<double>(fRinvnum); // Read matrix from file
        if( fabs(d - leg_nds[i])<tol){
	  RinvTemp = Rinvi[i];
	  break;
	}
        coef = w[i]/(d - leg_nds[i]);
	denom += coef;
	RinvTemp += coef*Rinvi[i];
      }
      RinvTemp = RinvTemp/denom;
    }
    else{
      BuildCompression(&Rtemp, &RinvTemp, x0, y0, x1, y1, radius, tol); // compute compression on-the-fly
    }
    
    const auto R_ = Rtemp;
    if(R != nullptr) (*R) = R_;
    const auto Rinv_ = RinvTemp;
    if(Rinv != nullptr) (*Rinv) = Rinv_;
  }
  

  template <class Real, Integer Order> void ICIP<Real,Order>::Setup() const {
    const auto& near_lst = disc_panels.GetNearList();
    const Long Nblocks = near_lst.Dim();
    if (icip_type_ != ICIPType::Adaptive && Nblocks != Kcorrec.Dim()) {
      Kcorrec.ReInit(Nblocks);
      Rprecon.ReInit(Nblocks);
      for (Long i = 0; i < Nblocks; i++) {
        Real x0, y0, x1, y1;
        std::tie(x0, y0) = disc_panels.DiscCoord(near_lst[i].disc_idx0);
        std::tie(x1, y1) = disc_panels.DiscCoord(near_lst[i].disc_idx1);

        if (icip_type_ == ICIPType::Compress) {
          GetPrecondBlock(nullptr, &Kcorrec[i], x0, y0, x1, y1, disc_panels.DiscRadius());
        } else if (icip_type_ == ICIPType::Precond) {
          GetPrecondBlock(&Rprecon[i], nullptr, x0, y0, x1, y1, disc_panels.DiscRadius());
        }
      }
    }
  }

  template <class Real, Integer Order> void ICIP<Real,Order>::ApplyMatrixBlocks(Vector<Real>& U, const Vector<Real>& F, const DiscPanelLst<Real,Order>& panel_lst, const Vector<typename DiscPanelLst<Real,Order>::NearData>& block_lst, const Vector<Matrix<Real>>& M_lst) {
    SCTL_ASSERT(block_lst.Dim() == M_lst.Dim());
    if (block_lst.Dim() == 0) return;

    const Long dof0 = F.Dim() / (panel_lst.Size()*Order);
    const Long dof1 = U.Dim() / (panel_lst.Size()*Order);
    SCTL_ASSERT(F.Dim() == panel_lst.Size()*Order*dof0);
    SCTL_ASSERT(U.Dim() == panel_lst.Size()*Order*dof1);
    U.SetZero();

    Matrix<Real> U_, F_; // temporary memory
    for (Long i = 0; i < block_lst.Dim(); i++) {
      const auto& block = block_lst[i];
      const Long offset_disc0 = (panel_lst.PanelIdxOffset(block.disc_idx0) + block.panel_idx_range0[0]) * Order;
      const Long offset_disc1 = (panel_lst.PanelIdxOffset(block.disc_idx1) + block.panel_idx_range1[0]) * Order;
      const Long N0 = (block.panel_idx_range0[1] - block.panel_idx_range0[0]) * Order;
      const Long N1 = (block.panel_idx_range1[1] - block.panel_idx_range1[0]) * Order;
      SCTL_ASSERT(M_lst[i].Dim(0) == (N0+N1)*dof0);
      SCTL_ASSERT(M_lst[i].Dim(1) == (N0+N1)*dof1);

      if (F_.Dim(0) != 1 || F_.Dim(1) != (N0+N1)*dof0) F_.ReInit(1, (N0+N1)*dof0);
      if (U_.Dim(0) != 1 || U_.Dim(1) != (N0+N1)*dof1) U_.ReInit(1, (N0+N1)*dof1);

      for (Long i = 0; i < N0*dof0; i++) F_[0][        i] = F[offset_disc0*dof0+i];
      for (Long i = 0; i < N1*dof0; i++) F_[0][N0*dof0+i] = F[offset_disc1*dof0+i];
      Matrix<Real>::GEMM(U_, F_, M_lst[i]);
      for (Long i = 0; i < N0*dof1; i++) U[offset_disc0*dof1+i] += U_[0][        i];
      for (Long i = 0; i < N1*dof1; i++) U[offset_disc1*dof1+i] += U_[0][N0*dof1+i];
    }
  }

  template <class Real, Integer Order> void ICIP<Real,Order>::ApplyBIOp(Vector<Real>* U, const Vector<Real>& sigma) const {
    const Long N = sigma.Dim();
    const auto& near_lst = this->disc_panels.GetNearList();
    if (icip_type_ == ICIPType::Precond && near_lst.Dim() > 0) {
      Setup();
      Vector<Real> R_sigma;
      ApplyPrecond(&R_sigma, sigma);
      this->ApplyBIOpDirect(U, R_sigma);
      { // U -= R_sigma_far
        Vector<Real> R_sigma_far, R_sigma_far_;
        Split(nullptr, &R_sigma_far, R_sigma);
        Merge(&R_sigma_far_, Vector<Real>(), R_sigma_far);
        (*U) -= R_sigma_far_;
      }
      (*U) += sigma;

    } else if (icip_type_ == ICIPType::Compress && near_lst.Dim() > 0) {
      Setup();
      Vector<Real> Ucorrec(N);
      this->ApplyMatrixBlocks(Ucorrec, sigma, this->disc_panels, this->disc_panels.GetNearList(), Kcorrec);
      this->ApplyBIOpDirect(U, sigma);
      (*U) += Ucorrec;

    } else {
      this->ApplyBIOpDirect(U, sigma);
    }
  }

  template <class Real, Integer Order> void ICIP<Real,Order>::ApplyPrecond(Vector<Real>* U, const Vector<Real>& sigma) const {
    const Long N = sigma.Dim();
    const auto& near_lst = this->disc_panels.GetNearList();
    if (icip_type_ == ICIPType::Precond && near_lst.Dim() > 0) {
      Setup();
      if (U->Dim() != N) U->ReInit(N);
      this->ApplyMatrixBlocks(*U, sigma, this->disc_panels, this->disc_panels.GetNearList(), Rprecon);
      { // U += sigma_far
        Vector<Real> sigma_far, sigma_far_;
        Split(nullptr, &sigma_far, sigma);
        Merge(&sigma_far_, Vector<Real>(), sigma_far);
        (*U) += sigma_far_;
      }

    } else {
      (*U) = sigma; // identity
    }
  }

  template <class Real, Integer Order> void ICIP<Real,Order>::SolveBIE(Vector<Real>& sigma, const Vector<Real>& rhs, const Real gmres_tol, const Long gmres_max_iter) const {
    const Long N = rhs.Dim();
    Vector<Real> rhs_ = rhs;
    SqrtScaling(rhs_);

    Vector<Real> sigma_;
    const auto BIOp = [this](Vector<Real>* U, const Vector<Real>& sigma) {
      auto sigma_ = sigma;
      InvSqrtScaling(sigma_);
      this->ApplyBIOp(U, sigma_);
      SqrtScaling(*U);
    };
    if (gmres_tol < 0) {
      Matrix<Real> M(N, N);
      { // Set M
        Vector<Real> x(N);
        for (Long i = 0; i < N; i++) {
          x.SetZero(); x[i] = 1;
          Vector<Real> U(N, M[i], false); U = 0;
          BIOp(&U, x);
          //std::cout<<i<<' '<<N<<'\n';
        }
        //std::string fname = std::string("M_") + (this->icip_type_==ICIPType::Adaptive ? "adaptive" : this->icip_type_==ICIPType::Compress ? "compress" : "precond") + ".mat";
        //if (std::is_same<Real,QuadReal>::value) M.Write(fname.c_str());
      }

      Matrix<Real> MM = Matrix<Real>(M).pinv(machine_eps<Real>());
      //std::string fname = std::string("Minv_") + (this->icip_type_==ICIPType::Adaptive ? "adaptive" : this->icip_type_==ICIPType::Compress ? "compress" : "precond") + ".mat";
      //if (std::is_same<Real,QuadReal>::value) MM.Write(fname.c_str());

      const Matrix<Real> rhs__(1, rhs_.Dim(), rhs_.begin());
      Matrix<Real> sigma = rhs__ * MM;
      sigma_ = Vector<Real>(sigma.Dim(1), sigma.begin(), false);
    } else {
      solver(&sigma_, BIOp, rhs_, gmres_tol, gmres_max_iter);
    }

    InvSqrtScaling(sigma_);
    this->ApplyPrecond(&sigma, sigma_);
  }

  template <class Real, Integer Order> void ICIP<Real,Order>::Split(Vector<Real>* v_near, Vector<Real>* v_far, const Vector<Real>& v) const {
    if (icip_type_ == ICIPType::Adaptive) {
      if (v_near && v_near->Dim()) v_near->ReInit(0);
      if (v_far) (*v_far) = v;
    } else {
      disc_panels.Split(v_near, v_far, v);
    }
  }

  template <class Real, Integer Order> void ICIP<Real,Order>::Merge(Vector<Real>* v, const Vector<Real>& v_near, const Vector<Real>& v_far) const {
    if (icip_type_ == ICIPType::Adaptive) {
      SCTL_ASSERT(!v_near.Dim());
      if (v) (*v) = v_far;
    } else {
      disc_panels.Merge(v, v_near, v_far);
    }
  }

  template <class Real, Integer Order> void ICIP<Real,Order>::SqrtScaling(Vector<Real>& v) const {
    const Long N = sqrt_wts.Dim();
    const Long dof = v.Dim() / N;
    SCTL_ASSERT(v.Dim() == N * dof);
    for (Long i = 0; i < N; i++) {
      const Real s = sqrt_wts[i];
      for (Long k = 0; k < dof; k++) v[i*dof+k] *= s;
    }
  }

  template <class Real, Integer Order> void ICIP<Real,Order>::InvSqrtScaling(Vector<Real>& v) const {
    const Long N = rsqrt_wts.Dim();
    const Long dof = v.Dim() / N;
    SCTL_ASSERT(v.Dim() == N * dof);
    for (Long i = 0; i < N; i++) {
      const Real s = rsqrt_wts[i];
      for (Long k = 0; k < dof; k++) v[i*dof+k] *= s;
    }
  }

  template <class Real, Integer Order> void ICIP<Real,Order>::disc_wise_outer_product(Vector<Real>& U, const Vector<Real>& sigma_far, const Vector<Real>& sigma_near, const Matrix<Real>& V0, const Matrix<Real>& V1) const {
    const auto add_interac = [this,&V0,&V1](Vector<Real>& U, const Vector<Real>& sigma, const Long disc_idx, const Long near_range0, const Long near_range1) {
      const Long Npanels = disc_panels.Size();
      const Long dof0 = V0.Dim(1) / Npanels;
      const Long dof1 = V1.Dim(1) / Npanels;
      const Long k0 = V0.Dim(0);
      SCTL_ASSERT(sigma.Dim() == Npanels * dof0);
      SCTL_ASSERT(V0.Dim(1) == Npanels * dof0);
      SCTL_ASSERT(V1.Dim(1) == Npanels * dof1);
      SCTL_ASSERT(V1.Dim(0) == k0);
      if (U.Dim() != Npanels * dof1) {
        U.ReInit(Npanels * dof1);
        U.SetZero();
      }

      Long range[2];
      const Long offset = disc_panels.PanelIdxOffset(disc_idx);
      const Long N = disc_panels.SurfWts(disc_idx).Dim() / Order;
      if (near_range1 > near_range0) {
        range[0] = offset + near_range0;
        range[1] = offset + near_range1;
      } else {
        range[0] = offset + 0;
        range[1] = offset + N;
      }

      for (Long k = 0; k < k0; k++) {
        Real sum = 0;
        for (Long j = range[0]*dof0; j < range[1]*dof0; j++) sum += V0[k][j] * sigma[j];

        if (near_range1 > near_range0) {
          for (Long j =   offset*dof1; j <   range[0]*dof1; j++) U[j] += sum * V1[k][j];
          for (Long j = range[1]*dof1; j < (offset+N)*dof1; j++) U[j] += sum * V1[k][j];
        } else {
          for (Long j =   offset*dof1; j < (offset+N)*dof1; j++) U[j] += sum * V1[k][j];
        }
      }
    };

    for (Long disc_idx = 0; disc_idx < disc_panels.DiscCount(); disc_idx++) {
      add_interac(U, sigma_far, disc_idx, 0, 0);
    }

    if (sigma_near.Dim()) {
      const auto& near_lst = disc_panels.GetNearList();
      for (const auto& near_block : near_lst) {
        add_interac(U, sigma_near, near_block.disc_idx0, near_block.panel_idx_range0[0], near_block.panel_idx_range0[1]);
        add_interac(U, sigma_near, near_block.disc_idx1, near_block.panel_idx_range1[0], near_block.panel_idx_range1[1]);
      }
    }
  }

}
