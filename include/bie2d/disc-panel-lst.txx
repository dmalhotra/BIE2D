namespace sctl {

  template <class Real, Integer Order> DiscPanelLst<Real,Order>::DiscPanelLst(const Comm& comm_) : comm(comm_){}

  template <class Real, Integer Order> DiscPanelLst<Real,Order>::~DiscPanelLst() {}

  template <class Real, Integer Order> void DiscPanelLst<Real,Order>::Init(const Vector<Real>& Xc_, const Real R_, bool adap, Real close_threshold) {
    SCTL_ASSERT(close_threshold < 2);
    Xc = Xc_;
    R = R_;

    const Long Ndisc = Xc.Dim() / COORD_DIM;
    theta_break.ReInit(Ndisc);
    for (auto& t : theta_break) t.ReInit(0);

    const Real pi = const_pi<Real>();
    const Real fac = 0.99;
    const Real dan = 2*asin<Real>(1/(2+close_threshold));
    const Real dhairya = acos<Real>(0.5+close_threshold/4);
    const Real close_angle = fac * std::min(dan, dhairya);
    const Real c = 0.5*close_angle;
    const Real max_size = 0.5*close_angle;

    Vector<Vector<Real>> breaks(Ndisc), endpoints(Ndisc);
    Vector<Vector<Long>> ireverse(Ndisc);
    Vector<Long> iforward;
    Long Npanel = 0;

    // Adaptively construct panel breakpoints on each disc
    // TODO: This is currently O(Ndisc^2), but should be O(Ndisc) with a kd-tree
    near_lst.ReInit(0);
    for (Long j = 0; j < Ndisc; j++) {
      for (Long k = j+1; k < Ndisc; k++) {

        const Real dx = Xc[2*j]   - Xc[2*k];
        const Real dy = Xc[2*j+1] - Xc[2*k+1];
        const Real dist = sqrt<Real>(dx*dx + dy*dy) - 2*R;
        SCTL_ASSERT(dist > 0);

        if (dist < close_threshold*R) {

          const Real theta_k = atan2<Real>(dy, dx); // in [-pi, pi]
          const Real theta_j = (theta_k < 0 ? theta_k+pi : theta_k-pi); // in [-pi, pi]

          // Add the endpoints of the near-region to the list of panel
          // breakpoints for each disc
          const Long start_idx_j = theta_break[j].Dim();
          const Long start_idx_k = theta_break[k].Dim();
          theta_break[j].PushBack(theta_j-c);
          theta_break[j].PushBack(theta_j+c);
          theta_break[k].PushBack(theta_k-c);
          theta_break[k].PushBack(theta_k+c);

          endpoints[j].PushBack(theta_j-c);
          endpoints[j].PushBack(theta_j+c);
          endpoints[k].PushBack(theta_k-c);
          endpoints[k].PushBack(theta_k+c);

          // Refine the near-region until the panel size is less than
          // the square root of the distance between the discs
          const Long nlevels = adap ? (Long)(2+0.5*log2<Real>(close_threshold*R/dist)) : 2;
          if (nlevels > 0) {
            theta_break[j].PushBack(theta_j);
            theta_break[k].PushBack(theta_k);
          }
          Real cl = c;
          for (Long l = 1; l < nlevels; l++) {
            cl = 0.5*cl;
            theta_break[j].PushBack(theta_j-cl);
            theta_break[j].PushBack(theta_j+cl);
            theta_break[k].PushBack(theta_k-cl);
            theta_break[k].PushBack(theta_k+cl);
          }

          NearData neardata;
          neardata.disc_idx0 = j;
          neardata.disc_idx1 = k;
          neardata.panel_idx_range0[0] = start_idx_j;
          neardata.panel_idx_range0[1] = start_idx_j+1; // the next breakpoint after start_idx is the near-region's end
          neardata.panel_idx_range1[0] = start_idx_k;
          neardata.panel_idx_range1[1] = start_idx_k+1;
          near_lst.PushBack(neardata);
        }
      }

      if (theta_break[j].Dim() == 0) {
        theta_break[j].PushBack(0);
        endpoints[j].PushBack(0);
        endpoints[j].PushBack(0);
      }

      // Now refine the leftover regions (if needed)
      std::sort(endpoints[j].begin(), endpoints[j].end());
      for (Long i = 1; i < endpoints[j].Dim(); i+=2) {
        const Real theta0 = endpoints[j][i];
        const Real theta1 = (i < endpoints[j].Dim()-1 ? endpoints[j][i+1] : endpoints[j][0]+2*pi);
        Real dt = theta1-theta0;
        const Long nuniform = (Long)ceil<Real>(dt/max_size);
        dt /= nuniform;
        for (Long l = 1; l < nuniform; l++) {
          theta_break[j].PushBack(theta0+l*dt);
        }
      }

      // Create index vector that puts the breaks in sorted order
      iforward.ReInit(theta_break[j].Dim());
      for (Long i = 0; i < iforward.Dim(); i++) {
        iforward[i] = i;
      }
      std::sort(iforward.begin(), iforward.end(), [&](Long ki, Long kj) {
        return theta_break[j][ki] < theta_break[j][kj];
      });

      // Create index vector that reverses the mapping
      ireverse[j].ReInit(theta_break[j].Dim());
      for (Long i = 0; i < theta_break[j].Dim(); i++) {
        ireverse[j][iforward[i]] = i;
      }

      // Do the actual sorting
      std::sort(theta_break[j].begin(), theta_break[j].end());

      // Add this disc's panels to the total number of panels
      Npanel += theta_break[j].Dim();
    }

    // Update the near interaction list with the sorted indices
    for (auto& neardata : near_lst) {
      neardata.panel_idx_range0[0] = ireverse[neardata.disc_idx0][neardata.panel_idx_range0[0]];
      neardata.panel_idx_range0[1] = ireverse[neardata.disc_idx0][neardata.panel_idx_range0[1]];
      neardata.panel_idx_range1[0] = ireverse[neardata.disc_idx1][neardata.panel_idx_range1[0]];
      neardata.panel_idx_range1[1] = ireverse[neardata.disc_idx1][neardata.panel_idx_range1[1]];
    }

    { // Init PanelLst, panel_cnt, panel_dsp
      Long offset = 0;
      panel_cnt.ReInit(Ndisc);
      panel_dsp.ReInit(Ndisc);
      Vector<Real> X(Npanel * Order * COORD_DIM);
      for (Long i = 0; i < Ndisc; i++) {
        const Long panel_count = theta_break[i].Dim();
        panel_cnt[i] = panel_count;
        panel_dsp[i] = offset;

        for (Long j = 0; j < panel_count; j++) {
          const Real theta0 = theta_break[i][j];
          const Real theta1 = (j < panel_count-1 ? theta_break[i][j+1] : theta_break[i][0]+2*const_pi<Real>());
          for (Long k = 0; k < Order; k++) {
            const Long node_idx = (offset + j) * Order + k;
            const Real theta = theta0 + (theta1-theta0) * PanelLstType::PanelNds()[k];
            X[node_idx*COORD_DIM+0] = Xc[i*COORD_DIM+0] + R * cos<Real>(theta);
            X[node_idx*COORD_DIM+1] = Xc[i*COORD_DIM+1] + R * sin<Real>(theta);
          }
        }
        offset += panel_count;
      }
      PanelLstType::Init(X);
    }

    { // Init X, Normal, SurfWts
      X.ReInit(Ndisc);
      Normal.ReInit(Ndisc);
      SurfWts__.ReInit(Ndisc);
      Long offset = 0;
      for (Long i = 0; i < Ndisc; i++) {
        const Long panel_count = theta_break[i].Dim();
        X[i].ReInit(panel_count*Order*COORD_DIM, (Iterator<Real>)PanelLstType::SurfCoord().begin() + offset*Order*COORD_DIM, false);
        Normal[i].ReInit(panel_count*Order*COORD_DIM, (Iterator<Real>)PanelLstType::SurfNormal().begin() + offset*Order*COORD_DIM, false);
        SurfWts__[i].ReInit(panel_count*Order, (Iterator<Real>)PanelLstType::SurfWts().begin() + offset*Order, false);
        offset += panel_count;
      }
    }

    { // setup split into near and far panels
      const auto& near_lst = GetNearList();
      const Long N = near_lst.Dim()*2;

      Vector<Long> panel_range(N*2+2);
      panel_range[0] = 0;
      panel_range.end()[-1] = this->Size();
      #pragma omp parallel for schedule(static)
      for (Long i = 0; i < near_lst.Dim(); i++) {
        const auto& n = near_lst[i];
        const Long panel_idx0 = PanelIdxOffset(n.disc_idx0);
        const Long panel_idx1 = PanelIdxOffset(n.disc_idx1);
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
  }

  template <class Real, Integer Order> Long DiscPanelLst<Real,Order>::DiscCount() const {
    return Xc.Dim() / COORD_DIM;
  }

  template <class Real, Integer Order> Real DiscPanelLst<Real,Order>::DiscRadius() const {
    return R;
  }

  template <class Real, Integer Order> std::tuple<Real,Real> DiscPanelLst<Real,Order>::DiscCoord(const Long disc_idx) const {
    return {Xc[disc_idx*COORD_DIM+0], Xc[disc_idx*COORD_DIM+1]};
  }

  template <class Real, Integer Order> const Vector<Real>& DiscPanelLst<Real,Order>::SurfCoord(const Long disc_idx) const {
    SCTL_ASSERT(disc_idx < X.Dim());
    return (disc_idx < 0 ? PanelLstType::SurfCoord() : X[disc_idx]);
  }
  template <class Real, Integer Order> const Vector<Real>& DiscPanelLst<Real,Order>::SurfNormal(const Long disc_idx) const {
    SCTL_ASSERT(disc_idx < Normal.Dim());
    return (disc_idx < 0 ? PanelLstType::SurfNormal() : Normal[disc_idx]);
  }
  template <class Real, Integer Order> const Vector<Real>& DiscPanelLst<Real,Order>::SurfWts(const Long disc_idx) const {
    SCTL_ASSERT(disc_idx < SurfWts__.Dim());
    return (disc_idx < 0 ? PanelLstType::SurfWts() : SurfWts__[disc_idx]);
  }
  template <class Real, Integer Order> const Vector<Real>& DiscPanelLst<Real,Order>::SurfTheta(const Long disc_idx) const {
    SCTL_ASSERT(disc_idx < theta_break.Dim());
    return (disc_idx < 0 ? theta_break_flat : theta_break[disc_idx]);
  }

  template <class Real, Integer Order> const Vector<typename DiscPanelLst<Real,Order>::NearData>& DiscPanelLst<Real,Order>::GetNearList() const {
    return near_lst;
  }

  template <class Real, Integer Order> Long DiscPanelLst<Real,Order>::PanelIdxOffset(const Long disc_idx) const {
    return panel_dsp[disc_idx];
  }

  template <class Real, Integer Order> void DiscPanelLst<Real,Order>::Split(Vector<Real>* v_near, Vector<Real>* v_far, const Vector<Real>& v) const {
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

  template <class Real, Integer Order> void DiscPanelLst<Real,Order>::Merge(Vector<Real>* v, const Vector<Real>& v_near, const Vector<Real>& v_far) const {
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

}
