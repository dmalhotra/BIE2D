#include <bie2d.hpp>

using namespace sctl;

template <class Real, class KerFn> void DoubleLayerTest(Real tol = 1e-10) {
  Vector<Disc<Real>> disc_lst(2);
  disc_lst[0] = Disc<Real>(0, 0, 0.5);
  disc_lst[1] = Disc<Real>(0, 1.5, 0.8);

  Vector<Real> X; // surface node coordinates
  GetGeom(disc_lst, &X);

  const Long N = NodeCount(disc_lst); // number of surface discretization nodes
  Vector<Real> F(N * KerFn::SrcDim());
  F = 1; // density

  Vector<Real> U;
  { // U <-- D[F]
    Matrix<Real> M, f_(1, F.Dim(), F.begin(), false);
    LayerPotentialMatrix<KerFn>(M, disc_lst, X, tol);
    if (U.Dim() != M.Dim(1)) U.ReInit(M.Dim(1));
    Matrix<Real> u_(1, U.Dim(), U.begin(), false);
    u_ = f_ * M;
  }

  auto max_norm = [](const sctl::Vector<Real>& U) {
    Real max_val = 0;
    for (const auto& u : U) max_val = std::max<Real>(max_val, sctl::fabs<Real>(u));
    return max_val;
  };
  std::cout<<"Error = "<< max_norm(U+0.5) <<'\n';
}

int main() {
  using Real = double;
  sctl::Profile::Enable(true);

  sctl::Profile::Tic("LaplaceDL");
  DoubleLayerTest<Real, Laplace2D_DxU>(1e-10);
  sctl::Profile::Toc();

  sctl::Profile::Tic("StokesDL");
  DoubleLayerTest<Real, Stokes2D_DxU>(1e-10);
  sctl::Profile::Toc();

  sctl::Profile::print();
  return 0;
}

