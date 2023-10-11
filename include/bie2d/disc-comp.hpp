#ifndef BIE2D_DISC_COMP_HPP
#define BIE2D_DISC_COMP_HPP

#include <sctl.hpp>
#include <bie2d/panel-lst.hpp>

template <class Real> class NearCorrection {

};


template <class Real, Integer Order, Integer digits> class DiscPanelLst : public PanelLst<Real,Order,digits,DiscPanelLst<Real,Order,digits>> {
  using PanelLstType = PanelLst<Real,Order,digits,Disc<Real,Order,digits>>;
  static constexpr Integer COORD_DIM = PanelLstType::COORD_DIM;

  public:

    struct NearData {
      Long disc_idx0, disc_idx1;
      StaticArray<Long,COORD_DIM> panel_idx_range0;
      StaticArray<Long,COORD_DIM> panel_idx_range1;
    };

    /**
     * Constructor
     */
    DiscPanelLst(const Comm& comm = Comm::Self());

    /**
     * @param[in] Xc coordinates of the center of each disc (length = Ndisc * COORD_DIM)
     *
     * @param[in] R radius of the discs.
     *
     * @param[in] adap whether to generate adaptive or coarse discretization.
     *
     * @param[in] far_dist_factor compute near interactions when distance is less than far_dist_factor*R.
     * The near panels makes an angle of min(acos(0.5+0.25*far_dist_factor), 2*asin(1/(2+far_dist_factor))).
     */
    void Init(const Vector<Real>& Xc_, const Real R, bool adap = false, Real far_dist_factor = 1.0) {
      Xc = Xc_;

      const Long Ndisc = Xc.Dim() / COORD_DIM;
      theta_break.ReInit(Ndisc);

      Long total_panel_count;
      theta_break_flat.ReInit(total_panel_count);
      Long offset = 0;
      for (Long i = 0; i < Ndisc; i++) {
        Long disc_i_panel_count;
        theta_break.ReInit(disc_i_panel_count, theta_break_flat.begin() + offset);
        offset += disc_i_panel_count;
      }
    }

    std::tuple<Real,Real> DiscCoord(const Long disc_idx) const {
      return {Xc[disc_idx*COORD_DIM+0], Xc[disc_idx*COORD_DIM+1]};
    }

    const Vector<Real>& SurfCoord(const Long disc_idx) const;
    const Vector<Real>& SurfNormal(const Long disc_idx) const;
    const Vector<Real>& SurfWts(const Long disc_idx) const;
    const Vector<Real>& SurfTheta(const Long disc_idx) const;
    //void GetGeom(Vector<Real>* X = nullptr, Vector<Real>* Normal = nullptr, Vector<Real>* SurfWts = nullptr) const;

    const Vector<NearData>& GetNearLst();

  private:

    Vector<Real> theta_break_flat;
    Vector<Vector<Real>> X, Normal, SurfWts, theta_break;

    Vector<NearData> near_lst;

    Vector<Real> Xc_;
};


template <class Real> class ICIP {
  public:

    ICIP();

    void Setup(const Vector<Real>& Xc, const Real R, bool icip = true);

  private:

    virtual void BuildPrecomp() = 0;
}

/**
 * Solve the Stokes mobility boundary integral equation.
 */
template <class Real> class DiscMobility {
  static constexpr Integer COORD_DIM = 2;
  public:

    DiscMobility(const Comm& comm = Comm::Self());

    //void Setup(const Vector<Real>& Xc, const Real R, bool icip = true);

    /**
     * @param[out] V velocity and angular velocity of each disc (length= = Nx3).
     *
     * @param[in] F force and torque on each disc (length = Nx3).
     *
     * @param[in] Vs slip velocity (length = Nnodes * COORD_DIM or zero slip if empty)
     */
    void SolveMobility(Vector<Real> V, const Vector<Real> F, const Vector<Real> Vs = Vector<Real>());

  private:

    virtual void BuildPrecomp(const Real R, const Long InterpOrder, const Real d_min, const Real d_max) {
      Vector<Real> Xc(2*COOR_DIM);
      Xc = 0;

      const auto& leg_nds = LegQuadRule<Real>::ComputeNds(InterpOrder);
      for (Long i = 0; i < InterpOrder; i++) { // loop over interpolation nodes
        const Real x = R + R/2 * d_min * exp(log(d_max/d_min) * leg_nds[i]);
        Xc[0*COORD_DIM+0] =-x;
        Xc[1*COORD_DIM+0] = x;
        DiscPanelLst<Real> panels_fine, panels_coarse;
        panels_coarse.Init(Xc, R, false, d_max);
        panels_fine.Init(Xc, R, true, d_max);

      }
    }

    Comm comm;
    DiscPanelLst<Real> panel_lst;
    NearCorrection<Real> near_correction; // stores the near quadrature corrections and the compressed near interactions
}



#endif
