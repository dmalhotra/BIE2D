#ifndef BIE2D_DISC_CAPACITANCE_HPP
#define BIE2D_DISC_CAPACITANCE_HPP

#include <sctl.hpp>
#include <bie2d/disc-icip.hpp>
#include <bie2d/kernels.hpp>

namespace sctl {

  /**
   * Solve the Laplace capacitance boundary integral equation.
   */
  template <class Real, Integer Order = 16> class DiscCapacitance : public ICIP<Real,Order> {
    static constexpr Integer COORD_DIM = 2;
    using ICIP_Base = ICIP<Real,Order>;

    public:

    DiscCapacitance(const Comm& comm_ = Comm::Self());

    virtual ~DiscCapacitance() = default;

    /**
     * Set all inputs before capacitance solve.
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
     * @param[out] Q charge on each disc (length= = N).
     *
     * @param[in] V potential at each conducting disc (length = N).
     *
     * @param[in] gmres_max_iter maximum number of GMRES iterations.
     */
    void Solve(Vector<Real>& Q, const Vector<Real> V, const Long gmres_max_iter);

    template <class RefReal> static void test(const Real R = 0.75, const Real eps = 1e-5, const Real tol = machine_eps<Real>(), const Long gmres_iter = 1000, const Comm comm = Comm::Self()) {
      const Long Ndisc = 2;
      Vector<Real> X(Ndisc*COORD_DIM), V(Ndisc);
      { // Set V, X
        X[0] = -(R+eps/2)*sqrt<Real>(0.5); X[1] = -(R+eps/2)*sqrt<Real>(0.5);
        X[2] =  (R+eps/2)*sqrt<Real>(0.5); X[3] =  (R+eps/2)*sqrt<Real>(0.5);

        V[0] = 0;
        V[1] = 1;
      }

      const auto real2quad = [](const Vector<Real>& v) {
        Vector<RefReal> w(v.Dim());
        for (Long i = 0; i < v.Dim(); i++) w[i] = (RefReal)v[i];
        return w;
      };

      Vector<RefReal> Q0; // reference solution
      DiscCapacitance<RefReal,Order+8> capacitance0(comm);
      capacitance0.Init(real2quad(X), R, tol*1e-2, ICIPType::Adaptive);
      capacitance0.Solve(Q0, real2quad(V), gmres_iter);

      Vector<Real> Q1, Q2, Q3;
      DiscCapacitance<Real,Order> capacitance(comm);
      capacitance.Init(X, R, tol, ICIPType::Adaptive);
      capacitance.Solve(Q1, V, gmres_iter);
      capacitance.Init(X, R, tol, ICIPType::Compress);
      capacitance.Solve(Q2, V, gmres_iter);
      capacitance.Init(X, R, tol, ICIPType::Precond );
      capacitance.Solve(Q3, V, gmres_iter);

      Real err_adap = 0, err_comp = 0, err_prec = 0, max_val = 0;
      for (const auto x : Q0) max_val = std::max<Real>(max_val, fabs((Real)x));
      for (const auto x : Q0-real2quad(Q1)) err_adap = std::max<Real>(err_adap, fabs((Real)x));
      for (const auto x : Q0-real2quad(Q2)) err_comp = std::max<Real>(err_comp, fabs((Real)x));
      for (const auto x : Q0-real2quad(Q3)) err_prec = std::max<Real>(err_prec, fabs((Real)x));
      std::cout<<"Error (adaptive) = "<<err_adap/max_val<<'\n';
      std::cout<<"Error (compress) = "<<err_comp/max_val<<'\n';
      std::cout<<"Error (precond)  = "<<err_prec/max_val<<'\n';
      std::cout<<Q0;
    }

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
     * Apply the Laplace capacitance operator directly on the current
     * discretization (direct -- without preconditioning).
     *
     * @param[out] U result vector.
     *
     * @param[in] F density vector.
     */
    virtual void ApplyBIOpDirect(Vector<Real>* U, const Vector<Real>& sigma) const;

    Laplace2D_FxU LaplaceSL_Ker; // Laplace SL kernel
    Laplace2D_DxU LaplaceDL_Ker; // Laplace DL kernel
    BoundaryIntegralOp<Real,Laplace2D_FxU> LaplaceSL_BIOp; // Laplace SL operator
    BoundaryIntegralOp<Real,Laplace2D_DxU> LaplaceDL_BIOp; // Laplace DL operator
    ParallelSolver<Real> solver; // GMRES solver
    Matrix<Real> V0; // rigid velocity basis for each disc (dimensions = 3 x Nnodes*COORD_DIM)
  };

}

#include <bie2d/disc-capacitance.txx>

#endif
