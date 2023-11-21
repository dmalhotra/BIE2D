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
    DiscPanelLst(const Comm& comm_ = Comm::Self());

    virtual ~DiscPanelLst();

    /**
     * Initialize.
     *
     * @param[in] Xc coordinates of the center of each disc (length = Ndisc * COORD_DIM)
     *
     * @param[in] R radius of the discs.
     *
     * @param[in] adap whether to generate adaptive or coarse discretization.
     *
     * @param[in] close_threshold compute near interactions when distance is less than close_threshold*R.
     * The near panels makes an angle of min(acos(0.5+0.25*close_threshold), 2*asin(1/(2+close_threshold))).
     */
    void Init(const Vector<Real>& Xc_, const Real R_, bool adap = false, Real close_threshold = 0.5);

    /**
     * Return the number of discs.
     */
    Long DiscCount() const;

    /**
     * Return the radius of the discs.
     */
    Real DiscRadius() const;

    /**
     * Get the coordinates of the center of a disc.
     *
     * @param[in] disc_idx index of the disc (0 <= disc_idx < DiscCount()).
     */
    std::tuple<Real,Real> DiscCoord(const Long disc_idx) const;

    /**
     * Get the surface discretization nodes.
     *
     * @param[in] disc_idx index of the disc (0 <= disc_idx < DiscCount()), or
     * if disc_idx=-1 then return the discretization nodes for all discs.
     *
     * @return const reference of the vector containing the surface
     * discretization nodes in AoS order (i.e. {x1,y1,x2,y2,...,xn,yn}).
     */
    const Vector<Real>& SurfCoord(const Long disc_idx) const;
    const Vector<Real>& SurfNormal(const Long disc_idx) const;
    const Vector<Real>& SurfWts(const Long disc_idx) const;
    const Vector<Real>& SurfTheta(const Long disc_idx) const;

    /**
     * @return const reference of a vector of NearData, containing the
     * near-interaction data.
     */
    const Vector<NearData>& GetNearList() const;

    /**
     * Returns the index of the first panel of a disc in the global panel list.
     *
     * @param[in] disc_idx index of the disc (0 <= disc_idx < DiscCount()).
     */
    Long PanelIdxOffset(const Long disc_idx) const;

    void Split(Vector<Real>* v_near, Vector<Real>* v_far, const Vector<Real>& v) const;

    void Merge(Vector<Real>* v, const Vector<Real>& v_near, const Vector<Real>& v_far) const;

    private:

    Comm comm;
    Real R;
    Vector<Real> Xc;

    Vector<Real> theta_break_flat;
    Vector<Vector<Real>> X, Normal, SurfWts__, theta_break;
    Vector<Long> panel_cnt, panel_dsp;

    Vector<NearData> near_lst;
    Vector<Long> near_dsp_orig, near_dsp, near_cnt;
    Vector<Long> far_dsp_orig, far_dsp, far_cnt;
  };

}

#include <bie2d/disc-panel-lst.txx>

#endif
