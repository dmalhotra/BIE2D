#ifndef BIE2D_DISC_COMP_HPP
#define BIE2D_DISC_COMP_HPP

template <class Real> class NearCorrection {

};

template <class Real, Integer Order, Integer digits> class DiscPanelLst : public PanelLst<Real,Order,digits,DiscPanelLst<Real,Order,digits>> {
  public:

    struct Neardata {
      StaticArray<Real,2> Xc0;
      StaticArray<Real,2> Xc1;
      StaticArray<Long,2> panel_idx_range0;
      StaticArray<Long,2> panel_idx_range1;
    };

    DiscPanelLst(const Comm& comm = Comm::Self());

    void Init(const Vector<Real>& X, const Real R, bool adap = false, Real far_dist_factor = 1.0) {
      Long total_panel_count, Ndisc;
      theta_break_flat.ReInit(total_panel_count);
      theta_break.ReInit(Ndisc);
      Long offset = 0;
      for (Long i = 0; i < Ndisc; i++) {
        Long disc_i_panel_count;
        theta_break.ReInit(disc_i_panel_count, theta_break_flat.begin() + offset);
        offset += disc_i_panel_count;
      }

      for (Long i = 0; i < Ndisc; i++) {
        for (Long j = 0; j < theta_break[i].Dim(); i++) {
          theta_break[i][j];
        }
      }

      for (const auto& theta_vec : theta_break) { // loop over discs
        for (const Real theta : theta_vec) {
          tehta;//
        }
      }

    }

    void GetGeom(Vector<Real>* X = nullptr, Vector<Real>* Normal = nullptr, Vector<Real>* SurfWts = nullptr) const;

    const Vector<NearData>& GetNearLst();

  private:

    Vector<Real> theta_break_flat;
    Vector<Vector<Real>> theta_break;

    Vector<NearData> near_lst;
};


template <class Real> class ICIP {
  public:

    ICIP();

    void Setup(const Vector<Real>& Xc, const Real R, bool icip = true);

  private:

    virtual void BuildPrecomp();
}

/**
 * Solve the Stokes mobility boundary integral equation.
 */
template <class Real> class DiscMobility {
  public:

    DiscMobility(const Comm& comm = Comm::Self());

    //void Setup(const Vector<Real>& Xc, const Real R, bool icip = true);

    /**
     * @param[out] V velocity and angular velocity of each disc (length= = Nx3).
     *
     * @param[in] F force and torque on each disc (length = Nx3).
     *
     * @param[in] Vs slip velocity (length = Nnodes * 2 or zero slip if empty)
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
