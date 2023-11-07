#ifndef BIE2D_DISC_MOBILITY_HPP
#define BIE2D_DISC_MOBILITY_HPP

#include <sctl.hpp>
#include <bie2d/disc-icip.hpp>
#include <bie2d/kernels.hpp>

namespace sctl {

  /**
   * Solve the Stokes mobility boundary integral equation.
   */
  template <class Real, Integer Order = 16> class DiscMobility : public ICIP<Real,Order> {
    static constexpr Integer COORD_DIM = 2;
    using ICIP_Base = ICIP<Real,Order>;

    public:

    DiscMobility(const Comm& comm_ = Comm::Self());

    virtual ~DiscMobility() = default;

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

    /**
     * @param[out] V velocity and angular velocity of each disc (length= = Nx3).
     *
     * @param[in] F force and torque on each disc (length = Nx3).
     *
     * @param[in] Vs slip velocity (length = Nnodes * COORD_DIM or zero slip if empty) (TODO: not implemented)
     *
     * @param[in] gmres_max_iter maximum number of GMRES iterations.
     */
    void Solve(Vector<Real>& V, const Vector<Real> F, const Vector<Real> Vs, const Long gmres_max_iter);

    private:

    /**
     * To be used for naming precomputed data files.
     */
    virtual const std::string& Name() const;

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
    virtual void BuildInteracBlock(Matrix<Real>& M, const DiscPanelLst<Real,Order> panel_lst, const typename DiscPanelLst<Real,Order>::NearData& interac_block, const Real tol) const;

    /**
     * Apply the Stokes mobility operator directly on the current
     * discretization (direct -- without preconditioning).
     *
     * @param[out] U result vector.
     *
     * @param[in] F density vector.
     */
    virtual void ApplyBIOpDirect(Vector<Real>* U, const Vector<Real>& sigma) const;

    /**
     * Construct the basis vectors for rigid body velocity of discs.
     */
    static void RigidVelocityBasis(Matrix<Real>& V, const DiscPanelLst<Real,Order>& disc_panels, const Long disc_idx = -1);

    Stokes2D_FxU StokesSL_Ker; // Stokes SL kernel
    Stokes2D_DxU StokesDL_Ker; // Stokes DL kernel
    BoundaryIntegralOp<Real,Stokes2D_FxU> StokesSL_BIOp; // Stokes SL operator
    BoundaryIntegralOp<Real,Stokes2D_FxU> StokesSL_BIOp_near, StokesSL_BIOp_far; // Stokes SL operator
    BoundaryIntegralOp<Real,Stokes2D_DxU> StokesDL_BIOp_near, StokesDL_BIOp_far; // Stokes DL operator
    ParallelSolver<Real> solver; // GMRES solver
    Matrix<Real> V0; // rigid velocity basis for each disc (dimensions = 3 x Nnodes*COORD_DIM)
  };

}

#include <bie2d/disc-mobility.txx>

#endif
