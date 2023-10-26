#ifndef BIE2D_DISC_PANEL_LST_HPP
#define BIE2D_DISC_PANEL_LST_HPP

#include <sctl.hpp>
#include <bie2d/panel-lst.hpp>

namespace sctl {

  /**
   * Derived class of PanelLst, for panelization of discs.
   */
  template <class Real, Integer Order = 16> class DiscPanelLst : public PanelLst<Real,Order> {
    using PanelLstType = PanelLst<Real,Order>;
    static constexpr Integer COORD_DIM = PanelLstType::CoordDim();

    public:

    struct NearData {
      Long disc_idx0, disc_idx1;
      StaticArray<Long,2> panel_idx_range0; // element range on disc0
      StaticArray<Long,2> panel_idx_range1; // element range on disc1
    };

    /**
     * Constructor
     */
    DiscPanelLst(const Comm& comm_ = Comm::Self()) : comm(comm_){}

    virtual ~DiscPanelLst() {}

    /**
     * Initialize.
     *
     * @param[in] Xc coordinates of the center of each disc (length = Ndisc * COORD_DIM)
     *
     * @param[in] R radius of the discs.
     *
     * @param[in] adap whether to generate adaptive or coarse discretization.
     *
     * @param[in] far_dist_factor compute near interactions when distance is less than far_dist_factor*R.
     * The near panels makes an angle of min(acos(0.5+0.25*far_dist_factor), 2*asin(1/(2+far_dist_factor))).
     */
    void Init(const Vector<Real>& Xc_, const Real R_, bool adap = false, Real far_dist_factor = 1.0) {
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

    /**
     * Return the number of discs.
     */
    Long DiscCount() const {
      return Xc.Dim() / COORD_DIM;
    }

    /**
     * Return the radius of the discs.
     */
    Real DiscRadius() const {
      return R;
    }

    /**
     * Get the coordinates of the center of a disc.
     *
     * @param[in] disc_idx index of the disc (0 <= disc_idx < DiscCount()).
     */
    std::tuple<Real,Real> DiscCoord(const Long disc_idx) const {
      return {Xc[disc_idx*COORD_DIM+0], Xc[disc_idx*COORD_DIM+1]};
    }

    /**
     * Get the surface discretization nodes.
     *
     * @param[in] disc_idx index of the disc (0 <= disc_idx < DiscCount()), or
     * if disc_idx=-1 then return the discretization nodes for all discs.
     *
     * @return const reference of the vector containing the surface
     * discretization nodes in AoS order (i.e. {x1,y1,x2,y2,...,xn,yn}).
     */
    const Vector<Real>& SurfCoord(const Long disc_idx) const {
      SCTL_ASSERT(disc_idx < X.Dim());
      return (disc_idx < 0 ? PanelLstType::SurfCoord() : X[disc_idx]);
    }
    const Vector<Real>& SurfNormal(const Long disc_idx) const {
      SCTL_ASSERT(disc_idx < Normal.Dim());
      return (disc_idx < 0 ? PanelLstType::SurfNormal() : Normal[disc_idx]);
    }
    const Vector<Real>& SurfWts(const Long disc_idx) const {
      SCTL_ASSERT(disc_idx < SurfWts__.Dim());
      return (disc_idx < 0 ? PanelLstType::SurfWts() : SurfWts__[disc_idx]);
    }
    const Vector<Real>& SurfTheta(const Long disc_idx) const {
      SCTL_ASSERT(disc_idx < theta_break.Dim());
      return (disc_idx < 0 ? theta_break_flat : theta_break[disc_idx]);
    }

    /**
     * @return const reference of a vector of NearData, containing the
     * near-interaction data.
     */
    const Vector<NearData>& GetNearList() const {
      return near_lst;
    }

    /**
     * Returns the index of the first panel of a disc in the global panel list.
     *
     * @param[in] disc_idx index of the disc (0 <= disc_idx < DiscCount()).
     */
    Long PanelIdxOffset(const Long disc_idx) const {
      return panel_dsp[disc_idx];
    }

    private:

    Comm comm;
    Real R;
    Vector<Real> Xc;

    Vector<Real> theta_break_flat;
    Vector<Vector<Real>> X, Normal, SurfWts__, theta_break;
    Vector<Long> panel_cnt, panel_dsp;

    Vector<NearData> near_lst;
  };

}

#endif
