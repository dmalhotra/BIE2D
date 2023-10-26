#ifndef BIE2D_KERNELS_HPP
#define BIE2D_KERNELS_HPP

#include <sctl.hpp>

namespace sctl {

  template <Integer digits, class Real, Integer N> Vec<Real,N> approx_log(const Vec<Real,N>& x) { // TODO: vectorize
    Vec<Real,N> logx;
    for (Integer i = 0; i < N; i++) logx.insert(i, (x[i] > 0 ? sctl::log<Real>(x[i]) : 0));
    return logx;
  }
  template <Integer digits, class Real, Integer N> Vec<Real,N> approx_inv(const Vec<Real,N>& x) { // TODO: vectorize
    Vec<Real,N> xrsqrt = approx_rsqrt<1>(x, x > Vec<Real,N>::Zero());
    return xrsqrt * xrsqrt;
  }

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
      u[0][0] = approx_log<digits>(r2);
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
      return 1 / (2 * const_pi<Real>());
    }
    template <Integer digits, class VecType> static void uKerMatrix(VecType (&u)[1][1], const VecType (&r)[2], const VecType (&n)[2], const void* ctx_ptr) {
      const VecType r2 = r[0]*r[0]+r[1]*r[1];
      const VecType r2inv = approx_inv<digits>(r2);
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
      const VecType r2inv = approx_inv<digits>(r2);
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
      using ScalarType = typename VecType::ScalarType;
      const VecType r2 = r[0]*r[0]+r[1]*r[1];
      const VecType r2inv = approx_inv<digits>(r2);
      const VecType log_rinv = approx_log<digits>(r2) * ((ScalarType)-0.5);
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
      return 1 / const_pi<Real>();
    }
    template <Integer digits, class VecType> static void uKerMatrix(VecType (&u)[2][2], const VecType (&r)[2], const VecType (&n)[2], const void* ctx_ptr) {
      const VecType r2 = r[0]*r[0]+r[1]*r[1];
      const VecType r4 = r2*r2;
      const VecType r4inv = approx_inv<digits>(r4);
      const VecType r_dot_n = r[0] * n[0] + r[1] * n[1];
      u[0][0] = r[0]*r[0]*r4inv*r_dot_n;
      u[0][1] = r[0]*r[1]*r4inv*r_dot_n;
      u[1][0] = r[1]*r[0]*r4inv*r_dot_n;
      u[1][1] = r[1]*r[1]*r4inv*r_dot_n;
    }
  };

  struct Stokes2D_FxT_ {
    static const std::string& Name() {
      static const std::string name = "Stokes2D-FxdU";
      return name;
    }
    static constexpr Integer FLOPS() {
      return 0;
    }
    template <class Real> static constexpr Real uKerScaleFactor() {
      return -1 / const_pi<Real>();
    }
    template <Integer digits, class VecType> static void uKerMatrix(VecType (&u)[2][4], const VecType (&r)[2], const VecType (&n)[2], const void* ctx_ptr) {
      const VecType r2 = r[0]*r[0]+r[1]*r[1];
      const VecType r2inv = approx_inv<digits>(r2);

      u[0][0] = r[0]*r[0]*r[0]*r2inv*r2inv;
      u[0][1] = r[0]*r[0]*r[1]*r2inv*r2inv;
      u[0][2] = r[0]*r[1]*r[0]*r2inv*r2inv;
      u[0][3] = r[0]*r[1]*r[1]*r2inv*r2inv;
      u[1][0] = r[1]*r[0]*r[0]*r2inv*r2inv;
      u[1][1] = r[1]*r[0]*r[1]*r2inv*r2inv;
      u[1][2] = r[1]*r[1]*r[0]*r2inv*r2inv;
      u[1][3] = r[1]*r[1]*r[1]*r2inv*r2inv;
    }
  };

  using Laplace2D_FxU = GenericKernel<Laplace2D_FxU_>;
  using Laplace2D_DxU = GenericKernel<Laplace2D_DxU_>;
  using Laplace2D_FxdU = GenericKernel<Laplace2D_FxdU_>;

  using Stokes2D_FxU = GenericKernel<Stokes2D_FxU_>;
  using Stokes2D_DxU = GenericKernel<Stokes2D_DxU_>;
  using Stokes2D_FxT = GenericKernel<Stokes2D_FxT_>;

}

#endif

