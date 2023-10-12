#ifndef BIE2D_DISC_COMP_HPP
#define BIE2D_DISC_COMP_HPP

#include <sctl.hpp>
#include <bie2d/panel-lst.hpp>

namespace sctl {

  template <class Real, Integer Order = 16> class DiscPanelLst : public PanelLst<Real,Order> {
    using PanelLstType = PanelLst<Real,Order>;
    static constexpr Integer COORD_DIM = PanelLstType::CoordDim();

    public:

    struct NearData {
      Long disc_idx0, disc_idx1;
      StaticArray<Long,COORD_DIM> panel_idx_range0;
      StaticArray<Long,COORD_DIM> panel_idx_range1;
    };

    /**
     * Constructor
     */
    DiscPanelLst(const Comm& comm_ = Comm::Self()) : comm(comm_){}

    virtual ~DiscPanelLst() {}

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
        R = 0.75;
        Xc.ReInit(0);
        Xc.PushBack(-0.8); Xc.PushBack(0);
        Xc.PushBack( 0.8); Xc.PushBack(0);

        theta_break_flat.ReInit(0);
        theta_break_flat.PushBack(-0.1/1);
        theta_break_flat.PushBack(-0.1/2);
        if (adap) {
          theta_break_flat.PushBack(-0.1/4);
          theta_break_flat.PushBack(-0.1/8);
          theta_break_flat.PushBack(-0.1/16);
          theta_break_flat.PushBack(-0.1/32);
          theta_break_flat.PushBack(-0.1/64);
        }
        theta_break_flat.PushBack(0);
        if (adap) {
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
        }
        theta_break_flat.PushBack( 0.5);
        if (adap) {
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
          Long disc_i_panel_count = (adap ? 22 : 12);
          theta_break[i].ReInit(disc_i_panel_count, theta_break_flat.begin() + offset, false);
          offset += disc_i_panel_count;
        }

        near_lst.ReInit(1);
        near_lst[0].disc_idx0 = 0;
        near_lst[0].disc_idx1 = 1;
        near_lst[0].panel_idx_range0[0] = 0;
        near_lst[0].panel_idx_range0[1] = 0 + (adap ? 14 : 4);
        near_lst[0].panel_idx_range0[0] = 4;
        near_lst[0].panel_idx_range0[1] = 4 + (adap ? 14 : 4);
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

    Long DiscCount() const {
      return Xc.Dim() / COORD_DIM;
    }

    Real DiscRadius() const {
      return R;
    }

    std::tuple<Real,Real> DiscCoord(const Long disc_idx) const {
      return {Xc[disc_idx*COORD_DIM+0], Xc[disc_idx*COORD_DIM+1]};
    }

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

    const Vector<NearData>& GetNearLst() const {
      return near_lst;
    }

    Long GetPanelIdx(const Long disc_idx, const Long disc_panel_idx) {
      return panel_dsp[disc_idx] + disc_panel_idx;
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

  //template <class Real, class Kernel> class DiscBIOp {
  //  static constexpr Integer COORD_DIM = Kernel.CoordDim();
  //  static constexpr Integer KDIM0 = Kernel.SrcDim();
  //  static constexpr Integer KDIM1 = Kernel.TrgDim();

  //  public:

  //  BoundaryIntegral(const Comm comm_) : comm(comm_) {}

  //  template <class PanelLstType> void Setup(const PanelLstType& panel_lst, const Vector<Real>& Xt = Vector<Real>()) {
  //    X = panel_lst.SurfCoord();
  //    Xn = panel_lst.SurfNormal();
  //    wts = panel_lst.SurfWts();
  //    Y = (Xt.Dim() ? Xt : X);
  //  }

  //  void EvaluatePotential(Vector<Real>& U, const Vector<Real>& F) const {
  //    const Long Nt = Y.Dim() / COORD_DIM;
  //    if (U.Dim() !=)
  //  }

  //  private:

  //  Comm comm;
  //  Vector<Real> X, Xn, wts; // src
  //  Vector<Real> Y; // trg
  //};


  template <class Real, Integer Order, class Derived> class ICIP {
    static constexpr Integer COORD_DIM = 2;

    public:

    ICIP(const Comm& comm_ = Comm::Self()) : comm(comm_) {}

    void Setup(const Vector<Real>& Xc, const Real R, bool icip = true) {
      disc_panels.Init(Xc, R, !icip);
    }

    private:

    virtual void BuildPrecomp() = 0;

    Comm comm;
    DiscPanelLst<Real,Order> disc_panels;
  };

  /**
   * Solve the Stokes mobility boundary integral equation.
   */
  template <class Real, Integer Order = 16> class DiscMobility : public ICIP<Real,Order,DiscMobility<Real,Order>> {
    using ICIPType = ICIP<Real,Order,DiscMobility<Real,Order>>;
    static constexpr Integer COORD_DIM = ICIPType::COORD_DIM;

    public:

    DiscMobility(const Comm& comm = Comm::Self());

    void Setup(const Vector<Real>& Xc, const Real R, bool icip = true) {
      ICIPType::Setup(Xc, R, icip);
    }

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
      Vector<Real> Xc(2*COORD_DIM);
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
    //BoundaryIntegral<Real,>
    //NearCorrection<Real> near_correction; // stores the near quadrature corrections and the compressed near interactions
  };


}

#endif
