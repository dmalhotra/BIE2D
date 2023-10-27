#ifndef BIE2D_DISC_ICIP_HPP
#define BIE2D_DISC_ICIP_HPP

#include <sctl.hpp>
#include <bie2d/disc-panel-lst.hpp>

namespace sctl {

  enum class ICIPType {
    Adaptive,
    Compress,
    Precond
  };

  /**
   * Base clas for ICIP.
   */
  template <class Real, Integer Order> class ICIP {
    static constexpr Integer InterpOrder = 64;
    static constexpr Integer COORD_DIM = 2;
    static constexpr Real d_max = 1.0; // largest distance between the discs is d_max*R
    static constexpr Real d_min = 1e-14; // smallest distance between the discs is d_min*R

    public:

    ICIP(const Comm& comm_ = Comm::Self());

    virtual ~ICIP();

    /**
     * Set all inputs before mobility solve.
     *
     * @param[in] Xc coordinates of the disc centers (length = Ndisc * COORD_DIM).
     *
     * @param[in] R radius of the discs.
     *
     * @param[in] tol accuracy tolerance.
     *
     * @param[in] icip_type adaptive, compressed or compress-precoditioned.
     */
    void Init(const Vector<Real>& Xc, const Real R, const Real tol, const ICIPType icip_type);

    protected:

    /**
     * To be used for naming precomputed data files.
     */
    virtual const std::string& Name() const = 0;

    /**
     * Build interaction block for a given set of panels.
     *
     * @param[out] M output interaction matrix block.
     *
     * @param[in] panel_lst panelization of the discs.
     *
     * @param[in] interac_block specifies the subset of panels of panel_lst for
     * which to construct the interaction block.
     *
     * @param[in] tol accuracy tolerance.
     */
    virtual void BuildInteracBlock(Matrix<Real>& M, const DiscPanelLst<Real,Order> panel_lst, const typename DiscPanelLst<Real,Order>::NearData& interac_block, const Real tol) const = 0;

    /**
     * Apply the boundary integral operator directly on the current
     * discretization (direct -- without preconditioning).
     *
     * @param[out] U result vector.
     *
     * @param[in] F density vector.
     */
    virtual void ApplyBIOpDirect(Vector<Real>* U, const Vector<Real>& F) const = 0;

    /**
     * Compute the compressed preconditioner matrix-block for interaction
     * between two discs of given radius with centers at (x0,y0) and (x1,y1).
     */
    void BuildCompression(Matrix<Real>* R, Matrix<Real>* Rinv, const Real x0, const Real y0, const Real x1, const Real y1, const Real radius, const Real tol) const;

    /**
     * Computed the compressed preconditioner from interpolation at
     * log-Legendre nodes. Cache matrices at interpolation nodes to files.
     */
    void GetPrecondBlock(Matrix<Real>* R, Matrix<Real>* Rinv, const Real x0, const Real y0, const Real x1, const Real y1, const Real radius, const Real tol) const;

    /**
     * Build the preconditioner and correction blocks (Kcorrec, Rprecon)
     */
    void Setup() const;

    /**
     * Apply a sparse matrix that is represented by a set of small matrix-blocks.
     *
     * @param[out] U result vector, it must be preallocated to the appropriate length.
     *
     * @param[in] F density vector.
     *
     * @param[in] panel_lst geometry description.
     *
     * @param[in] block_lst vector of NearData, describin the matrix sparsity pattern.
     *
     * @param[in] M_lst list of small matrix-blocks.
     */
    static void ApplyMatrixBlocks(Vector<Real>& U, const Vector<Real>& F, const DiscPanelLst<Real,Order> panel_lst, const Vector<typename DiscPanelLst<Real,Order>::NearData>& block_lst, const Vector<Matrix<Real>>& M_lst);

    /**
     * Apply the boundary integral operator.
     */
    void ApplyBIOp(Vector<Real>* U, const Vector<Real>& sigma) const;

    /**
     * Apply the block diagonal preconditioner.
     */
    void ApplyPrecond(Vector<Real>* U, const Vector<Real>& sigma) const;

    Comm comm;
    Real tol_;
    ICIPType icip_type_;
    DiscPanelLst<Real,Order> disc_panels;
    mutable Vector<Matrix<Real>> Kcorrec; // blocks to add to Kc
    mutable Vector<Matrix<Real>> Rprecon; // block diagonal precond
  };

  template <class Real, Integer Order> constexpr Real ICIP<Real,Order>::d_max;
  template <class Real, Integer Order> constexpr Real ICIP<Real,Order>::d_min;

}

#include <bie2d/disc-icip.txx>

#endif
