#include <bie2d.hpp>

using namespace sctl;

template <class Real> Real max_norm(const sctl::Vector<Real>& U) {
  Real max_val = 0;
  for (const auto& u : U) max_val = std::max<Real>(max_val, sctl::fabs<Real>(u));
  return max_val;
};

template <class Real, class KerSL, class KerDL, class KerGrad> void GreensIdentityTest(Real tol = 1e-10) {
  const Real R = 1.0;
  const Long Npanels = 20;
  Vector<Disc<Real>> disc_lst(1);
  disc_lst[0] = Disc<Real>(0.0, 0.0, R, Npanels);

  Vector<Real> X, Xn; // surface node coordinates and normals
  GetGeom(disc_lst, &X, &Xn);

  Vector<Real> X0(2), F0(KerSL::SrcDim()); // source point coordinates and charge
  X0[0] = -0.9/sqrt(2.0);
  X0[1] = -0.9/sqrt(2.0);
  F0 = 1;

  const Long N = NodeCount(disc_lst);
  Vector<Real> U0, dU0, dU0_dn(N*KerSL::SrcDim()); // boundary data
  KerSL::KerEval(U0, X, X0, Vector<Real>(), F0);
  KerGrad::KerEval(dU0, X, X0, Vector<Real>(), F0);
  for (Long i = 0; i < N; i++) {
    for (Long k = 0; k < KerSL::TrgDim(); k++) {
      dU0_dn[i*KerSL::SrcDim()+k] = dU0[(i*KerSL::TrgDim()+k)*2+0] * Xn[i*2+0] + dU0[(i*KerSL::TrgDim()+k)*2+1] * Xn[i*2+1];
    }
  }

  { // Green's Identity (off-surface)
    const Long Nt = Npanels*16;
    Vector<Real> Xt(Nt * 2); // exterior target points
    for (Long i = 0; i < Nt; i++) {
      const Real x[2] = {drand48()-0.5, drand48()-0.5};
      const Real r = sqrt<Real>(x[0]*x[0] + x[1]*x[1]);
      Xt[i*2+0] = x[0] / r * (drand48()*0.01 + R);
      Xt[i*2+1] = x[1] / r * (drand48()*0.01 + R);
    }

    Vector<Real> Uref; // reference potential at Xt
    KerSL::KerEval(Uref, Xt, X0, Vector<Real>(), F0);

    Vector<Real> Usl, Udl;
    LayerPotential<KerSL>(Usl, disc_lst, Xt, dU0_dn, tol);
    LayerPotential<KerDL>(Udl, disc_lst, Xt, U0, tol);
    Vector<Real> U = (Udl - Usl);
    std::cout<<"Error (off-surface) = "<< max_norm(U - Uref) / max_norm(Uref) <<'\n';
  }

  { // Green's Identity (on-surface)
    Vector<Real> Usl, Udl;
    LayerPotential<KerSL>(Usl, disc_lst, X, dU0_dn, tol);
    LayerPotential<KerDL>(Udl, disc_lst, X, U0, tol);
    Vector<Real> U = (U0*0.5 + Udl) - Usl;
    std::cout<<"Error (on-surface) = "<< max_norm(U - U0) / max_norm(U0) <<'\n';
  }
}

int main() {
  using Real = long double;
  sctl::Profile::Enable(true);

  sctl::Profile::Tic("GreensIdentity");
  GreensIdentityTest<Real, Laplace2D_FxU, Laplace2D_DxU, Laplace2D_FxdU>();
  GreensIdentityTest<Real, Stokes2D_FxU, Stokes2D_DxU, Stokes2D_FxT>();
  sctl::Profile::Toc();

  sctl::Profile::print();
  return 0;
}

