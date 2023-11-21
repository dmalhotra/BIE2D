#ifndef BIE2D_DISC_BIOP_HPP
#define BIE2D_DISC_BIOP_HPP

#include <sctl.hpp>
#include <bie2d/disc-panel-lst.hpp>

namespace sctl {

  template <class Real, Integer Order, class Kernel> class DiscBIOp {
    static constexpr Integer COORD_DIM   = Kernel::CoordDim();
    static constexpr Integer KDIM0 = Kernel::SrcDim();
    static constexpr Integer KDIM1 = Kernel::TrgDim();

    public:

      DiscBIOp(const Comm& comm_);

      void Init(const DiscPanelLst<Real,Order>& disc_panels_, bool exclude_near_, const Real tol_);

      void ComputePotential(Vector<Real>& U, const Vector<Real>& F) const;

    private:

      void Split(Vector<Real>* v_near, Vector<Real>* v_far, const Vector<Real>& v) const;

      void Merge(Vector<Real>* v, const Vector<Real>& v_near, const Vector<Real>& v_far) const;

      static void ApplyMatrixBlocks(Vector<Real>& U, const Vector<Real>& F, const DiscPanelLst<Real,Order>& panel_lst, const Vector<typename DiscPanelLst<Real,Order>::NearData>& block_lst, const Vector<Matrix<Real>>& M_lst);

      Comm comm;
      Kernel ker;

      DiscPanelLst<Real,Order> disc_panels;
      bool exclude_near;
      Real tol;

      Vector<typename DiscPanelLst<Real,Order>::NearData> near_lst;
      Vector<Long> near_dsp_orig, near_dsp, near_cnt;
      Vector<Long> far_dsp_orig, far_dsp, far_cnt;

      Vector<Real> Xnear, Xfar;
      PanelLst<Real,Order> panels_near, panels_far;
      BoundaryIntegralOp<Real,Kernel> biop_n2f, biop_f2a;

      mutable Vector<Matrix<Real>> K_near;
  };

}

#include <bie2d/disc-biop.txx>

#endif
