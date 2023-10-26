namespace sctl {

  template <class Real, Integer Order> DiscPanelLst<Real,Order>::DiscPanelLst(const Comm& comm_) : comm(comm_){}

  template <class Real, Integer Order> DiscPanelLst<Real,Order>::~DiscPanelLst() {}

  template <class Real, Integer Order> void DiscPanelLst<Real,Order>::Init(const Vector<Real>& Xc_, const Real R_, bool adap, Real far_dist_factor) {
    SCTL_ASSERT(far_dist_factor < 2);
    Xc = Xc_;
    R = R_;

    const Long Ndisc = Xc.Dim() / COORD_DIM;
    theta_break.ReInit(Ndisc);
    { // Set theta_break_flat, theta_break
      Long total_panel_count = 0;
      theta_break_flat.ReInit(total_panel_count);
      Long offset = 0;
      for (Long i = 0; i < Ndisc; i++) {
        Long disc_i_panel_count = 0;
        theta_break[i].ReInit(disc_i_panel_count, theta_break_flat.begin() + offset, false);
        offset += disc_i_panel_count;
      }
    }

    if (1) { // init 2-disc case for testing
      if (R <= 0) R = 0.75;
      if (Xc.Dim() != 2*COORD_DIM) {
        Xc.ReInit(0);
        Xc.PushBack(-0.8); Xc.PushBack(0);
        Xc.PushBack( 0.8); Xc.PushBack(0);
      }

      theta_break_flat.ReInit(0);
      theta_break_flat.PushBack(-0.1/1);
      theta_break_flat.PushBack(-0.1/2);
      if (adap) {
        theta_break_flat.PushBack(-0.1/4);
        theta_break_flat.PushBack(-0.1/8);
        theta_break_flat.PushBack(-0.1/16);
        theta_break_flat.PushBack(-0.1/32);
        theta_break_flat.PushBack(-0.1/64);
        theta_break_flat.PushBack(-0.1/128);
        theta_break_flat.PushBack(-0.1/256);
        theta_break_flat.PushBack(-0.1/512);
        theta_break_flat.PushBack(-0.1/1024);
        theta_break_flat.PushBack(-0.1/2048);
        theta_break_flat.PushBack(-0.1/4096);
      }
      theta_break_flat.PushBack(0);
      if (adap) {
        theta_break_flat.PushBack( 0.1/4096);
        theta_break_flat.PushBack( 0.1/2048);
        theta_break_flat.PushBack( 0.1/1024);
        theta_break_flat.PushBack( 0.1/512);
        theta_break_flat.PushBack( 0.1/256);
        theta_break_flat.PushBack( 0.1/128);
        theta_break_flat.PushBack( 0.1/64);
        theta_break_flat.PushBack( 0.1/32);
        theta_break_flat.PushBack( 0.1/16);
        theta_break_flat.PushBack( 0.1/8);
        theta_break_flat.PushBack( 0.1/4);
      }
      theta_break_flat.PushBack( 0.1/2);
      theta_break_flat.PushBack( 0.1/1);
      theta_break_flat.PushBack( 0.2);
      theta_break_flat.PushBack( 0.3);
      theta_break_flat.PushBack( 0.4);
      theta_break_flat.PushBack( 0.5);
      theta_break_flat.PushBack( 0.6);
      theta_break_flat.PushBack( 0.7);
      theta_break_flat.PushBack( 0.8);

      theta_break_flat.PushBack( 0.0);
      theta_break_flat.PushBack( 0.1);
      theta_break_flat.PushBack( 0.2);
      theta_break_flat.PushBack( 0.3);
      theta_break_flat.PushBack( 0.5 - 0.1/1);
      theta_break_flat.PushBack( 0.5 - 0.1/2);
      if (adap) {
        theta_break_flat.PushBack( 0.5 - 0.1/4);
        theta_break_flat.PushBack( 0.5 - 0.1/8);
        theta_break_flat.PushBack( 0.5 - 0.1/16);
        theta_break_flat.PushBack( 0.5 - 0.1/32);
        theta_break_flat.PushBack( 0.5 - 0.1/64);
        theta_break_flat.PushBack( 0.5 - 0.1/128);
        theta_break_flat.PushBack( 0.5 - 0.1/256);
        theta_break_flat.PushBack( 0.5 - 0.1/512);
        theta_break_flat.PushBack( 0.5 - 0.1/1024);
        theta_break_flat.PushBack( 0.5 - 0.1/2048);
        theta_break_flat.PushBack( 0.5 - 0.1/4096);
      }
      theta_break_flat.PushBack( 0.5);
      if (adap) {
        theta_break_flat.PushBack( 0.5 + 0.1/4096);
        theta_break_flat.PushBack( 0.5 + 0.1/2048);
        theta_break_flat.PushBack( 0.5 + 0.1/1024);
        theta_break_flat.PushBack( 0.5 + 0.1/512);
        theta_break_flat.PushBack( 0.5 + 0.1/256);
        theta_break_flat.PushBack( 0.5 + 0.1/128);
        theta_break_flat.PushBack( 0.5 + 0.1/64);
        theta_break_flat.PushBack( 0.5 + 0.1/32);
        theta_break_flat.PushBack( 0.5 + 0.1/16);
        theta_break_flat.PushBack( 0.5 + 0.1/8);
        theta_break_flat.PushBack( 0.5 + 0.1/4);
      }
      theta_break_flat.PushBack( 0.5 + 0.1/2);
      theta_break_flat.PushBack( 0.5 + 0.1/1);
      theta_break_flat.PushBack( 0.7);
      theta_break_flat.PushBack( 0.8);
      theta_break_flat.PushBack( 0.9);
      theta_break_flat *= 2*const_pi<Real>();

      const Long Ndisc = Xc.Dim() / COORD_DIM;
      theta_break.ReInit(Ndisc);
      Long offset = 0;
      for (Long i = 0; i < Ndisc; i++) {
        Long disc_i_panel_count = (adap ? 34 : 12);
        theta_break[i].ReInit(disc_i_panel_count, theta_break_flat.begin() + offset, false);
        offset += disc_i_panel_count;
      }

      near_lst.ReInit(1);
      near_lst[0].disc_idx0 = 0;
      near_lst[0].disc_idx1 = 1;
      near_lst[0].panel_idx_range0[0] = 0;
      near_lst[0].panel_idx_range0[1] = 0 + (adap ? 26 : 4);
      near_lst[0].panel_idx_range1[0] = 4;
      near_lst[0].panel_idx_range1[1] = 4 + (adap ? 26 : 4);
    }

    { // Init PanelLst, panel_cnt, panel_dsp
      const Long Ndisc = Xc.Dim() / COORD_DIM;
      const Long Npanel = theta_break_flat.Dim();

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

}
