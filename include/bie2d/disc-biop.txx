namespace sctl {

  template <class Real, Integer Order, class Kernel> DiscBIOp<Real,Order,Kernel>::DiscBIOp(const Comm& comm_) : comm(comm_), biop_n2f(ker, false, comm), biop_f2a(ker, false, comm) {}

  template <class Real, Integer Order, class Kernel> void DiscBIOp<Real,Order,Kernel>::Init(const DiscPanelLst<Real,Order>& disc_panels_, bool exclude_near_, const Real tol_) {
    disc_panels = disc_panels_;
    exclude_near = exclude_near_;
    tol = tol_;

    // Near-panels are the inner-panels (excluding two side panels on
    // each disc) of all close-to-toching blocks and far-panels includes everything else.
    { // Setup split into near and far panels.
      near_lst = disc_panels.GetNearList();
      for (auto& n : near_lst) { // exclude the two side panels from each disc
        n.panel_idx_range0[0]++;
        n.panel_idx_range0[1]--;
        n.panel_idx_range1[0]++;
        n.panel_idx_range1[1]--;
      }
      const Long N = near_lst.Dim()*2;

      Vector<Long> panel_range(N*2+2);
      panel_range[0] = 0;
      panel_range.end()[-1] = disc_panels.Size();
      #pragma omp parallel for schedule(static)
      for (Long i = 0; i < near_lst.Dim(); i++) {
        const auto& n = near_lst[i];
        const Long panel_idx0 = disc_panels.PanelIdxOffset(n.disc_idx0);
        const Long panel_idx1 = disc_panels.PanelIdxOffset(n.disc_idx1);
        panel_range[1+i*4+0] = panel_idx0 + n.panel_idx_range0[0];
        panel_range[1+i*4+1] = panel_idx0 + n.panel_idx_range0[1];
        panel_range[1+i*4+2] = panel_idx1 + n.panel_idx_range1[0];
        panel_range[1+i*4+3] = panel_idx1 + n.panel_idx_range1[1];
      }
      omp_par::merge_sort(panel_range.begin()+1, panel_range.end()-1);

      far_dsp_orig.ReInit(N+1);
      far_dsp.ReInit(N+1);
      far_cnt.ReInit(N+1);
      far_dsp[0] = 0;

      near_dsp_orig.ReInit(N);
      near_dsp.ReInit(N);
      near_cnt.ReInit(N);
      near_dsp[0] = 0;

      for (Long i = 0; i < N; i++) {
        const Long n0 = panel_range[i*2+0];
        const Long n1 = panel_range[i*2+1];
        const Long n2 = panel_range[i*2+2];
        far_dsp_orig[i] = n0;
        far_cnt[i] = n1-n0;
        near_dsp_orig[i] = n1;
        near_cnt[i] = n2-n1;
      }
      far_dsp_orig[N] = panel_range[N*2+0];
      far_cnt[N] = panel_range[N*2+1] - panel_range[N*2+0];
      omp_par::scan(far_cnt.begin(), far_dsp.begin(), N+1);
      omp_par::scan(near_cnt.begin(), near_dsp.begin(), N);
    }
    K_near.ReInit(0);

    const auto& X = disc_panels.SurfCoord(-1);
    { // Set Xnear, Xfar, panels_near, panels_far
      if (exclude_near) {
        Split(&Xnear, &Xfar, X);
      } else {
        Xnear.ReInit(0);
        Xfar = X;
      }
      panels_near.Init(Xnear);
      panels_far.Init(Xfar);
    }

    biop_n2f.AddElemList(panels_near);
    biop_n2f.SetTargetCoord(Xfar);
    biop_n2f.SetAccuracy(tol);

    biop_f2a.AddElemList(panels_far);
    biop_f2a.SetTargetCoord(X);
    biop_f2a.SetAccuracy(tol);
  }

  template <class Real, Integer Order, class Kernel> void DiscBIOp<Real,Order,Kernel>::ComputePotential(Vector<Real>& U, const Vector<Real>& F) const {
    const Long Npanels = disc_panels.Size();
    SCTL_ASSERT(F.Dim() == Npanels * Order * KDIM0);
    if (U.Dim() != Npanels * Order * KDIM1) U.ReInit(Npanels * Order * KDIM1);
    U.SetZero();

    if (exclude_near) {
      Vector<Real> F_near, F_far, U_, U__;
      Split(&F_near, &F_far, F);
      biop_f2a.ComputePotential(U, F_far);
      biop_n2f.ComputePotential(U_, F_near);
      Merge(&U__, Vector<Real>(), U_);
      U += U__;

      { // Subtract contributions to/from the two side panels of each close-to-touching block.
        const auto& near_lst = disc_panels.GetNearList();
        const Long Nblocks = near_lst.Dim();
        if (K_near.Dim() != Nblocks) {
          K_near.ReInit(Nblocks);
          for (Long i = 0; i < Nblocks; i++) {
            const Long disc_offset0 = disc_panels.PanelIdxOffset(near_lst[i].disc_idx0);
            const Long disc_offset1 = disc_panels.PanelIdxOffset(near_lst[i].disc_idx1);
            const Long start0 = disc_offset0 + near_lst[i].panel_idx_range0[0];
            const Long   end0 = disc_offset0 + near_lst[i].panel_idx_range0[1];
            const Long start1 = disc_offset1 + near_lst[i].panel_idx_range1[0];
            const Long   end1 = disc_offset1 + near_lst[i].panel_idx_range1[1];
            const Long N0 = end0-start0;
            const Long N1 = end1-start1;

            const auto Xt = disc_panels.SurfCoord(-1);
            Vector<Real> Xt_((N0+N1)*Order*COORD_DIM);
            for (Long j = 0; j < N0*Order*COORD_DIM; j++) Xt_[                   j] = Xt[start0*Order*COORD_DIM+j];
            for (Long j = 0; j < N1*Order*COORD_DIM; j++) Xt_[N0*Order*COORD_DIM+j] = Xt[start1*Order*COORD_DIM+j];

            K_near[i].ReInit((N0+N1)*Order*KDIM0, (N0+N1)*Order*KDIM1);
            K_near[i].SetZero();
            Matrix<Real> K0(N0*Order*KDIM0, (N0+N1)*Order*KDIM1, K_near[i].begin(), false);
            Matrix<Real> K1(N1*Order*KDIM0, (N0+N1)*Order*KDIM1, K_near[i].begin() + (N0*Order*KDIM0)*((N0+N1)*Order*KDIM1), false);
            disc_panels.template LayerPotentialMatrix<Kernel,-1>(K0, Xt_, tol, start0, end0);
            disc_panels.template LayerPotentialMatrix<Kernel,-1>(K1, Xt_, tol, start1, end1);

            #pragma omp parallel for schedule(static)
            for (Long j = Order*KDIM0; j < (N0-1)*Order*KDIM0; j++) {
              for (Long k = Order*KDIM1; k < (N0-1)*Order*KDIM1; k++) {
                K0[j][k] = 0;
              }
              for (Long k = (N0+1)*Order*KDIM1; k < (N0+N1-1)*Order*KDIM1; k++) {
                K0[j][k] = 0;
              }
            }
            #pragma omp parallel for schedule(static)
            for (Long j = Order*KDIM0; j < (N1-1)*Order*KDIM0; j++) {
              for (Long k = Order*KDIM1; k < (N0-1)*Order*KDIM1; k++) {
                K1[j][k] = 0;
              }
              for (Long k = (N0+1)*Order*KDIM1; k < (N0+N1-1)*Order*KDIM1; k++) {
                K1[j][k] = 0;
              }
            }
          }
        }

        Vector<Real> U_near(U.Dim()); U_near.SetZero();
        ApplyMatrixBlocks(U_near, F, disc_panels, near_lst, K_near);
        U -= U_near;
      }

      { // Add ineractions between near panels of different close-to-touching blocks (using direct GL-quadrature)
        const auto& X = disc_panels.SurfCoord(-1);
        const auto& Xn = disc_panels.SurfNormal(-1);
        const auto& W = disc_panels.SurfWts(-1);

        Vector<Real> F_(Npanels * Order * KDIM0);
        for (Long i = 0; i < Npanels; i++) {
          for (Long j = 0; j < Order; j++) {
            const Real w = W[i*Order+j];
            for (Long k = 0; k < KDIM0; k++) {
              F_[(i*Order+j)*KDIM0+k] = F[(i*Order+j)*KDIM0+k] * w;
            }
          }
        }

        Vector<Real> U0, U1, U_; // temporary memory
        for (Long i = 0; i < near_lst.Dim(); i++) {
          const auto& near_block = near_lst[i];
          const Long offset_disc0 = (disc_panels.PanelIdxOffset(near_block.disc_idx0) + near_block.panel_idx_range0[0]) * Order;
          const Long offset_disc1 = (disc_panels.PanelIdxOffset(near_block.disc_idx1) + near_block.panel_idx_range1[0]) * Order;
          const Long N0 = (near_block.panel_idx_range0[1] - near_block.panel_idx_range0[0]) * Order;
          const Long N1 = (near_block.panel_idx_range1[1] - near_block.panel_idx_range1[0]) * Order;

          const Vector<Real> Xs0(N0*COORD_DIM, (Iterator<Real>)X.begin() + offset_disc0*COORD_DIM, false);
          const Vector<Real> Xs1(N1*COORD_DIM, (Iterator<Real>)X.begin() + offset_disc1*COORD_DIM, false);
          const Vector<Real> Xn0(N0*COORD_DIM, (Iterator<Real>)Xn.begin() + offset_disc0*COORD_DIM, false);
          const Vector<Real> Xn1(N1*COORD_DIM, (Iterator<Real>)Xn.begin() + offset_disc1*COORD_DIM, false);
          const Vector<Real> F0_(N0*KDIM0, F_.begin() + offset_disc0*KDIM0, false);
          const Vector<Real> F1_(N1*KDIM0, F_.begin() + offset_disc1*KDIM0, false);

          U0.SetZero();
          U1.SetZero();
          Kernel::template Eval<Real,true>(U0, Xnear, Xs0, Xn0, F0_, -1, Ptr2ConstItr<char>(&ker,sizeof(ker)));
          Kernel::template Eval<Real,true>(U1, Xnear, Xs1, Xn1, F1_, -1, Ptr2ConstItr<char>(&ker,sizeof(ker)));
          U0 += U1;

          Merge(&U_, U0, Vector<Real>());
          Vector<Real> U0_(N0*KDIM1, U_.begin() + offset_disc0*KDIM1, false);
          Vector<Real> U1_(N1*KDIM1, U_.begin() + offset_disc1*KDIM1, false);
          U0_.SetZero();
          U1_.SetZero();
          U += U_;
        }
      }
    } else {
      biop_f2a.ComputePotential(U, F);
    }
  }

  template <class Real, Integer Order, class Kernel> void DiscBIOp<Real,Order,Kernel>::Split(Vector<Real>* v_near, Vector<Real>* v_far, const Vector<Real>& v) const {
    if ((v_near == nullptr && v_far == nullptr) || !far_cnt.Dim()) return;
    const Long N = far_dsp_orig.end()[-1] +  far_cnt.end()[-1];
    const Long Nfar  =  far_dsp.end()[-1] +  far_cnt.end()[-1];
    const Long Nnear = near_dsp.end()[-1] + near_cnt.end()[-1];
    const Long dof = v.Dim() / N;
    SCTL_ASSERT(v.Dim() == N*dof);

    if (v_near) {
      if (v_near->Dim() != Nnear*dof) v_near->ReInit(Nnear*dof);
      #pragma omp parallel for schedule(static)
      for (Long i = 0; i < near_cnt.Dim(); i++) {
        const Long offset_orig = near_dsp_orig[i];
        const Long offset = near_dsp[i];
        for (Long j = 0; j < near_cnt[i]; j++) {
          for (Long k = 0; k < dof; k++) {
            (*v_near)[(offset+j)*dof+k] = v[(offset_orig+j)*dof+k];
          }
        }
      }
    }
    if (v_far) {
      if (v_far->Dim() != Nfar*dof) v_far->ReInit(Nfar*dof);
      #pragma omp parallel for schedule(static)
      for (Long i = 0; i < far_cnt.Dim(); i++) {
        const Long offset_orig = far_dsp_orig[i];
        const Long offset = far_dsp[i];
        for (Long j = 0; j < far_cnt[i]; j++) {
          for (Long k = 0; k < dof; k++) {
            (*v_far)[(offset+j)*dof+k] = v[(offset_orig+j)*dof+k];
          }
        }
      }
    }
  }

  template <class Real, Integer Order, class Kernel> void DiscBIOp<Real,Order,Kernel>::Merge(Vector<Real>* v, const Vector<Real>& v_near, const Vector<Real>& v_far) const {
    if (v == nullptr || !far_cnt.Dim()) return;
    const Long N = far_dsp_orig.end()[-1] +  far_cnt.end()[-1];
    const Long Nfar  =  far_dsp.end()[-1] +  far_cnt.end()[-1];
    const Long Nnear = near_dsp.end()[-1] + near_cnt.end()[-1];
    const Long dof = (v_near.Dim() ? v_near.Dim()/Nnear : v_far.Dim()/Nfar);
    if (v_near.Dim()) SCTL_ASSERT(v_near.Dim() == Nnear*dof);
    if ( v_far.Dim()) SCTL_ASSERT( v_far.Dim() ==  Nfar*dof);

    if (v->Dim() != N*dof) v->ReInit(N*dof);
    if (v_near.Dim()) {
      #pragma omp parallel for schedule(static)
      for (Long i = 0; i < near_cnt.Dim(); i++) {
        const Long offset_orig = near_dsp_orig[i];
        const Long offset = near_dsp[i];
        for (Long j = 0; j < near_cnt[i]; j++) {
          for (Long k = 0; k < dof; k++) {
            (*v)[(offset_orig+j)*dof+k] = v_near[(offset+j)*dof+k];
          }
        }
      }
    } else {
      #pragma omp parallel for schedule(static)
      for (Long i = 0; i < near_cnt.Dim(); i++) {
        const Long offset_orig = near_dsp_orig[i]*dof;
        for (Long j = 0; j < near_cnt[i]*dof; j++) {
          (*v)[offset_orig+j] = 0;
        }
      }
    }
    if (v_far .Dim()) {
      #pragma omp parallel for schedule(static)
      for (Long i = 0; i < far_cnt.Dim(); i++) {
        const Long offset_orig = far_dsp_orig[i];
        const Long offset = far_dsp[i];
        for (Long j = 0; j < far_cnt[i]; j++) {
          for (Long k = 0; k < dof; k++) {
            (*v)[(offset_orig+j)*dof+k] = v_far[(offset+j)*dof+k];
          }
        }
      }
    } else {
      #pragma omp parallel for schedule(static)
      for (Long i = 0; i < far_cnt.Dim(); i++) {
        const Long offset_orig = far_dsp_orig[i]*dof;
        for (Long j = 0; j < far_cnt[i]*dof; j++) {
          (*v)[offset_orig+j] = 0;
        }
      }
    }
  }

  template <class Real, Integer Order, class Kernel> void DiscBIOp<Real,Order,Kernel>::ApplyMatrixBlocks(Vector<Real>& U, const Vector<Real>& F, const DiscPanelLst<Real,Order>& panel_lst, const Vector<typename DiscPanelLst<Real,Order>::NearData>& block_lst, const Vector<Matrix<Real>>& M_lst) {
    SCTL_ASSERT(block_lst.Dim() == M_lst.Dim());
    if (block_lst.Dim() == 0) return;

    const Long dof0 = F.Dim() / (panel_lst.Size()*Order);
    const Long dof1 = U.Dim() / (panel_lst.Size()*Order);
    SCTL_ASSERT(F.Dim() == panel_lst.Size()*Order*dof0);
    SCTL_ASSERT(U.Dim() == panel_lst.Size()*Order*dof1);

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

}
