#ifndef BIE2D_KERNELS_HPP
#define BIE2D_KERNELS_HPP

#include <sctl.hpp>

namespace sctl {

  template <class uKernel> class GenericKernel2D : public uKernel {
    template <class VecType, Integer K0, Integer K1, Integer D, class ...T> static constexpr Integer get_DIM  (void (*uKer)(VecType (&u)[K0][K1], const VecType (&r)[D], T... args)) { return D;  }
    template <class VecType, Integer K0, Integer K1, Integer D, class ...T> static constexpr Integer get_KDIM0(void (*uKer)(VecType (&u)[K0][K1], const VecType (&r)[D], T... args)) { return K0; }
    template <class VecType, Integer K0, Integer K1, Integer D, class ...T> static constexpr Integer get_KDIM1(void (*uKer)(VecType (&u)[K0][K1], const VecType (&r)[D], T... args)) { return K1; }

    static constexpr Integer DIM   = get_DIM  (uKernel::template uKerMatrix<0,double>);
    static constexpr Integer KDIM0 = get_KDIM0(uKernel::template uKerMatrix<0,double>);
    static constexpr Integer KDIM1 = get_KDIM1(uKernel::template uKerMatrix<0,double>);

    public:

    static constexpr Integer CoordDim() {
      return DIM;
    }
    static constexpr Integer SrcDim() {
      return KDIM0;
    }
    static constexpr Integer TrgDim() {
      return KDIM1;
    }

    private:

    friend uKernel;
  };

  struct Laplace2D_FxU_ {
    static const std::string& Name() {
      static const std::string name = "Laplace2D-FxU";
      return name;
    }
    static constexpr Integer FLOPS() {
      return 0;
    }
    template <class Real> static constexpr Real uKerScaleFactor() {
      return -1 / (4 * const_pi<Real>());
    }
    template <Integer digits, class VecType> static void uKerMatrix(VecType (&u)[1][1], const VecType (&r)[2], const VecType (&n)[2], const void* ctx_ptr) {
      const VecType r2 = r[0]*r[0]+r[1]*r[1];
      u[0][0] = (r2 > 0 ? sctl::log<VecType>(r2) : (VecType)0);
    }
  };

  struct Laplace2D_DxU_ {
    static const std::string& Name() {
      static const std::string name = "Laplace2D-DxU";
      return name;
    }
    static constexpr Integer FLOPS() {
      return 0;
    }
    template <class Real> static constexpr Real uKerScaleFactor() {
      return -1 / (2 * const_pi<Real>());
    }
    template <Integer digits, class VecType> static void uKerMatrix(VecType (&u)[1][1], const VecType (&r)[2], const VecType (&n)[2], const void* ctx_ptr) {
      const VecType r2 = r[0]*r[0]+r[1]*r[1];
      const VecType r2inv = (r2 > 0 ? 1/r2 : (VecType)0);
      const VecType r_dot_n = r[0] * n[0] + r[1] * n[1];
      u[0][0] = r2inv * r_dot_n;
    }
  };

  struct Laplace2D_FxdU_ {
    static const std::string& Name() {
      static const std::string name = "Laplace2D-FxdU";
      return name;
    }
    static constexpr Integer FLOPS() {
      return 0;
    }
    template <class Real> static constexpr Real uKerScaleFactor() {
      return -1 / (2 * const_pi<Real>());
    }
    template <Integer digits, class VecType> static void uKerMatrix(VecType (&u)[1][2], const VecType (&r)[2], const VecType (&n)[2], const void* ctx_ptr) {
      const VecType r2 = r[0]*r[0]+r[1]*r[1];
      const VecType r2inv = (r2 > 0 ? 1/r2 : (VecType)0);
      u[0][0] = r[0] * r2inv;
      u[0][1] = r[1] * r2inv;
    }
  };

  struct Stokes2D_FxU_ {
    static const std::string& Name() {
      static const std::string name = "Stokes2D-FxU";
      return name;
    }
    static constexpr Integer FLOPS() {
      return 0;
    }
    template <class Real> static constexpr Real uKerScaleFactor() {
      return 1 / (4 * const_pi<Real>());
    }
    template <Integer digits, class VecType> static void uKerMatrix(VecType (&u)[2][2], const VecType (&r)[2], const VecType (&n)[2], const void* ctx_ptr) {
      const VecType r2 = r[0]*r[0]+r[1]*r[1];
      const VecType r2inv = (r2 > 0 ? 1/r2 : 0);
      const VecType log_rinv = 0.5*sctl::log(r2inv);
      u[0][0] = log_rinv + r[0]*r[0]*r2inv;
      u[0][1] =            r[0]*r[1]*r2inv;
      u[1][0] =            r[1]*r[0]*r2inv;
      u[1][1] = log_rinv + r[1]*r[1]*r2inv;
    }
  };

  struct Stokes2D_DxU_ {
    static const std::string& Name() {
      static const std::string name = "Stokes2D-DxU";
      return name;
    }
    static constexpr Integer FLOPS() {
      return 0;
    }
    template <class Real> static constexpr Real uKerScaleFactor() {
      return -1 / const_pi<Real>();
    }
    template <Integer digits, class VecType> static void uKerMatrix(VecType (&u)[2][2], const VecType (&r)[2], const VecType (&n)[2], const void* ctx_ptr) {
      const VecType r2 = r[0]*r[0]+r[1]*r[1];
      const VecType r4inv = (r2 > 0 ? 1/(r2*r2) : 0);
      const VecType r_dot_n = r[0] * n[0] + r[1] * n[1];
      u[0][0] = r[0]*r[0]*r4inv*r_dot_n;
      u[0][1] = r[0]*r[1]*r4inv*r_dot_n;
      u[1][0] = r[1]*r[0]*r4inv*r_dot_n;
      u[1][1] = r[1]*r[1]*r4inv*r_dot_n;
    }
  };

  // TODO
  //struct Stokes2D_FxdU_ {
  //  static const std::string& Name() {
  //    static const std::string name = "Stokes2D-FxdU";
  //    return name;
  //  }
  //  static constexpr Integer FLOPS() {
  //    return 0;
  //  }
  //  template <class Real> static constexpr Real uKerScaleFactor() {
  //    static_assert(false, "Not implemented");
  //    return 0;
  //  }
  //  template <Integer digits, class VecType> static void uKerMatrix(VecType (&u)[2][4], const VecType (&r)[2], const VecType (&n)[2], const void* ctx_ptr) {
  //    static_assert(false, "Not implemented");
  //  }
  //};

  using Laplace2D_FxU = GenericKernel2D<Laplace2D_FxU_>;
  using Laplace2D_DxU = GenericKernel2D<Laplace2D_DxU_>;
  using Laplace2D_FxdU = GenericKernel2D<Laplace2D_FxdU_>;

  using Stokes2D_FxU = GenericKernel2D<Stokes2D_FxU_>;
  using Stokes2D_DxU = GenericKernel2D<Stokes2D_DxU_>;
  //using Stokes2D_FxdU = GenericKernel2D<Stokes2D_FxdU_>; // TODO

}

#endif
