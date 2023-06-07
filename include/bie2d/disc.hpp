#ifndef BIE2D_DISC_HPP
#define BIE2D_DISC_HPP

#include <sctl.hpp>
#include <bie2d/panel-lst.hpp>

namespace sctl {

  /**
   * Implements layer potentials on discs. The discs are discretized using
   * Gauss-Legendre panels, where each panel has "Order" nodes.
   */
  template <class Real, Integer Order = 16, Integer digits = 10> class Disc : private PanelLst<Real,Order,digits,Disc<Real,Order,digits>> {
    using PanelType = PanelLst<Real,Order,digits,Disc<Real,Order,digits>>;
    static constexpr Integer COORD_DIM = PanelType::COORD_DIM;

    public:

    static constexpr Integer CoordDim();

    /**
     * Returns vector of length Order with Gauss-Legendre nodes for the interval (0,1).
     */
    static const Vector<Real>& PanelNds();

    /**
     * Returns vector of length Order with Gauss-Legendre weights for the interval (0,1).
     */
    static const Vector<Real>& PanelWts();


    Disc();

    /**
     * \brief Initialize a Disc object with Nunif uniform panels.
     */
    Disc(const Real x, const Real y, const Real radius, const long Nunif = 8);

    /**
     * \brief Initialize a Disc object with panels specified by starting and
     * ending theta values in theta0 and theta1 respectively.
     */
    Disc(const Real x, const Real y, const Real radius, const Vector<Real>& theta0, const Vector<Real>& theta1);

    /**
     * \brief Initialize a Disc object with panels specified by starting and
     * ending theta values in theta0 and theta1 respectively.
     */
    void Init(const Real x, const Real y, const Real radius, const Vector<Real>& theta0, const Vector<Real>& theta1);



    Real Radius() const;

    Real Coord(int i) const;

    /**
     * \brief Return the number of panels.
     */
    Long PanelCount() const;

    /**
     * \brief Return the number of discretization nodes
     * = PanelCount() * Order
     */
    Long NodeCount() const;

    /**
     * \brief Return the starting and ending theta values of a panel.
     */
    std::pair<Real,Real> PanelRange(Long panel_idx) const;

    /**
     * \brief Return the position vector X, the normal vector Normal, and the
     * theta values at the surface discretization nodes. The values are in AoS
     * order.
     */
    void GetGeom(Vector<Real>* X = nullptr, Vector<Real>* Normal = nullptr, Vector<Real>* theta = nullptr) const;



    /**
     * Compute the surface integral from function values F given at surface
     * discretization nodes.
     */
    void BoundaryIntegralDirect(Vector<Real>& I, const Vector<Real>& F) const;

    /**
     * Evaluate the potential at the target points Xt, resulting from the
     * convolution of a kernel function with a surface density F given at the
     * surface discretization nodes.
     */
    template <class KerFn> void LayerPotential(Vector<Real>& U, const Vector<Real>& Xt, const Vector<Real>& F, const Real tol = 1e-10) const;

    /**
     * Return the operator matrix M for the layer potential operator.
     * M has the dimensions Nt*K0 x Ns*K1
     * where Nt is the number of target points, Ns is the number of surface
     * discretization points, and K0 x K1 is the kernel function dimension (eg.
     * Laplace single-layer is 1x1 and Laplace gradient is 1x2).
     */
    template <class KerFn> void LayerPotentialMatrix(Matrix<Real>& M, const Vector<Real>& Xt, const Real tol = 1e-10) const;

    private:

    Real radius_;
    StaticArray<Real,COORD_DIM> coord_;
    Vector<Real> theta0_, theta1_;
  };

  // Operations on Vector<Disc>

  /**
   * \brief Return the total number of discretization nodes for a vector of discs.
   */
  template <class Disc> Long NodeCount(const Vector<Disc>& disc_lst);

  /**
   * \brief Return the position vector X, the normal vector Normal, and the
   * theta values at the surface discretization nodes for all discs in disc_lst.
   * The values are in AoS order.
   */
  template <class Real, class Disc> void GetGeom(const Vector<Disc>& disc_lst, Vector<Real>* X = nullptr, Vector<Real>* Normal = nullptr);

  /**
   * Compute the surface integral from function values F given at surface
   * discretization nodes on all discs in disc_lst.
   */
  template <class Real, class Disc> void BoundaryIntegralDirect(Vector<Real>& I, const Vector<Disc>& disc_lst, const Vector<Real>& F);

  /**
   * Evaluate the potential at the target points Xt, resulting from the
   * convolution of a kernel function with a surface density F given at the
   * surface discretization nodes of all discs in disc_lst.
   */
  template <class Real, class KerFn, class Disc> void LayerPotentialMatrix(Matrix<Real>& M, const Vector<Disc>& disc_lst, const Vector<Real>& Xt, const Real tol = 1e-10);

  /**
   * Return the operator matrix M for the layer potential operator.
   * M has the dimensions Nt*K0 x Ns*K1
   * where Nt is the number of target points, Ns is the total number of surface
   * discretization points on all discs in disc_lst, and K0 x K1 is the kernel
   * function dimension.
   */
  template <class Real, class KerFn, class Disc> void LayerPotential(Vector<Real>& U, const Vector<Disc>& disc_lst, const Vector<Real>& Xt, const Vector<Real>& F, const Real tol = 1e-10);

}

#include <bie2d/disc.txx>

#endif
