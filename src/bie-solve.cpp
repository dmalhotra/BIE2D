#include <bie2d.hpp>

using namespace sctl;

template <class Real> void SolveCapacitance(const Vector<Disc<Real>>& disc_lst, const Real tol = 1e-10, const Long max_iter = 200) {
  using KerSL = Laplace2D_FxU;
  using KerDL = Laplace2D_DxU;

  Vector<Real> Xt;
  GetGeom(disc_lst, &Xt);

  // Setup boundary integral operator: K = 1/2 + D + S (exterior limit)
  auto BIOp = [&disc_lst,&Xt,&tol](Vector<Real>* Ax, const Vector<Real>& x) {
    Vector<Real> Usl, Udl;
    {
      Matrix<Real> Msl, Mdl;
      LayerPotentialMatrix<KerSL>(Msl, disc_lst, Xt, tol);
      LayerPotentialMatrix<KerDL>(Mdl, disc_lst, Xt, tol);

      Usl.ReInit(Msl.Dim(1));
      Udl.ReInit(Mdl.Dim(1));
      Matrix<Real> Usl_(1, Usl.Dim(), Usl.begin(), false);
      Matrix<Real> Udl_(1, Udl.Dim(), Udl.begin(), false);
      const Matrix<Real> x_(1, x.Dim(), (Iterator<Real>)x.begin(), false);
      Usl_ = x_ * Msl;
      Udl_ = x_ * Mdl;
    }
    (*Ax) = 0.5*x + Udl + Usl;
  };

  Vector<Real> U0; // Dirichlet boundary conditions
  for (Long i = 0; i < disc_lst.Dim(); i++) { // Set U0
    const Long N = disc_lst[i].NodeCount();
    for (Long j = 0; j < N; j++) {
      U0.PushBack(i); // U0 = i for i-th disc
    }
  }

  Vector<Real> sigma; // unknown density
  ParallelSolver<Real> solver; // Linear solver (GMRES)
  solver(&sigma, BIOp, U0, tol, max_iter);

  // Print total charge on each disc
  Long offset = 0;
  for (auto& disc : disc_lst) {
    const Long N = disc.NodeCount();
    Vector<Real> sigma_(N, sigma.begin() + offset, false);

    Vector<Real> q;
    disc.BoundaryIntegralDirect(q, sigma_);
    std::cout<<q;

    offset += N;
  }
}

int main() {
  using Real = double;
  sctl::Profile::Enable(true);

  Vector<Disc<Real>> disc_lst(2);
  disc_lst[0] = Disc<Real>(0, 0, 0.499, 6);
  disc_lst[1] = Disc<Real>(0, 1, 0.499, 6);

  sctl::Profile::Tic("Capacitance");
  SolveCapacitance<>(disc_lst);
  sctl::Profile::Toc();

  sctl::Profile::print();
  return 0;
}

