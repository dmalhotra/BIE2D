namespace sctl {

  template <class Real, Integer Order> DiscPanelLst<Real,Order>::DiscPanelLst(const Comm& comm_) : comm(comm_){}

  template <class Real, Integer Order> DiscPanelLst<Real,Order>::~DiscPanelLst() {}

  template <class Real, Integer Order> void DiscPanelLst<Real,Order>::Init(const Vector<Real>& Xc_, const Real R_, bool adap, Real far_dist_factor) {
    SCTL_ASSERT(far_dist_factor < 2);
    Xc = Xc_;
    R = R_;

    const Long Ndisc = Xc.Dim() / COORD_DIM;
    theta_break.ReInit(Ndisc);

    const Real pi = const_pi<Real>();
    Real fac = 1;
    Real close_threshold = 0.5;
    Real dan = 2*asin(1/(2+close_threshold));
    Real dhairya = acos(0.5+close_threshold/4);
    Real close_angle = fac * std::min(dan, dhairya);
    Real c = 0.5*close_angle;

    Vector<Vector<Real>> breaks(Ndisc);
    Vector<Vector<Long>> ireverse(Ndisc);
    Vector<Long> iforward;
    Long Npanel = 0;

    // Adaptively construct panel breakpoints on each disc
    // TODO: This is currently O(Ndisc^2), but should be O(Ndisc) with a kd-tree
    for (Long j = 0; j < Ndisc; j++)
    {
        for (Long k = j+1; k < Ndisc; k++)
        {
            Real dx = Xc[2*j]   - Xc[2*k];
            Real dy = Xc[2*j+1] - Xc[2*k+1];
            Real dist = sqrt(dx*dx + dy*dy) - 2*R;
            if ( dist < close_threshold*R )
            {
                Real theta_j = fmod(atan2(dy,dx)+pi, 2*pi);
                Real theta_k = fmod(theta_j+pi, 2*pi);

                // Add the endpoints of the near-region to the list of panel
                // breakpoints for each disc
                Long start_idx_j = theta_break[j].Dim();
                Long start_idx_k = theta_break[k].Dim();
                theta_break[j].PushBack(theta_j-c);
                theta_break[j].PushBack(theta_j+c);
                theta_break[k].PushBack(theta_k-c);
                theta_break[k].PushBack(theta_k+c);

                // Refine the near-region until the panel size is less than
                // the square root of the distance between the discs
                Long nlevels = adap ? (Long)ceil(log2(close_angle/sqrt(dist))) : 4;
                if (nlevels > 0)
                {
                    theta_break[j].PushBack(theta_j);
                    theta_break[k].PushBack(theta_k);
                }
                Real cl = c;
                for (Long l = 1; l < nlevels; l++)
                {
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
                neardata.panel_idx_range0[1] = theta_break[j].Dim()-1;
                neardata.panel_idx_range1[0] = start_idx_k;
                neardata.panel_idx_range1[1] = theta_break[k].Dim()-1;
                near_lst.PushBack(neardata);
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
      neardata.panel_idx_range0[0] = ireverse[neardata.disc_idx0][neardata.panel_idx_range0[0]] * Order;
      neardata.panel_idx_range0[1] = ireverse[neardata.disc_idx0][neardata.panel_idx_range0[1]] * Order;
      neardata.panel_idx_range1[0] = ireverse[neardata.disc_idx1][neardata.panel_idx_range1[0]] * Order;
      neardata.panel_idx_range1[1] = ireverse[neardata.disc_idx1][neardata.panel_idx_range1[1]] * Order;
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
