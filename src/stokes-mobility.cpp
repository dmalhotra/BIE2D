#include <bie2d.hpp>
using namespace sctl;

const char data_path[]="data2/";
constexpr Integer COORD_DIM = 2;

template <Integer ElemOrder, class Real> Vector<Real> DiscMobilitySolve(const Vector<Real>& F, const Vector<Real>& X, const Real R, const ICIPType icip_type, const Real tol = 1e-10) {
  const Comm comm = Comm::Self();
  DiscMobility<Real,ElemOrder> disc_mobility(comm);

  Vector<Real> V;
  disc_mobility.Init(X, R, tol, icip_type);
  disc_mobility.Solve(V, F, Vector<Real>(), 300);
  return V;
}

int main() {
  using Real = double;
  const Real tol = machine_eps<Real>()*64;

  const Real dt = 0.1;
  const Real R = 0.75;
  const Long Ndisc = 2;

  Vector<Real> X(Ndisc*COORD_DIM), F(Ndisc*3);
  { // Set F, X
    X = 0;
    X[0] = -0.7501; X[1] = 0.0;
    X[2] =  0.7501; X[3] = 0.0;
    //X[4] =  0.0; X[5] =-0.5;

    F = 0;
    F[0] = 0; F[1] = 1; F[2] = 0;
    F[3] = 0; F[4] = 0; F[5] = 0;
    //F[6] = 0; F[7] = 0; F[8] = 0;
  }

  //for (Long i = 0; i < 1000; i++) {
    const Vector<Real> V0 = DiscMobilitySolve<16>(F, X, R, ICIPType::Adaptive, tol);
    const Vector<Real> V1 = DiscMobilitySolve<16>(F, X, R, ICIPType::Compress, tol);
    const Vector<Real> V2 = DiscMobilitySolve<16>(F, X, R, ICIPType::Precond, tol);
    Real err_compress = 0, err_precond = 0, max_val = 0;
    for (const auto x : V0) max_val = std::max<Real>(max_val, fabs(x));
    for (const auto x : V0-V1) err_compress = std::max<Real>(err_compress, fabs(x));
    for (const auto x : V0-V2) err_precond = std::max<Real>(err_precond, fabs(x));
    std::cout<<"Error (compress) = "<<err_compress/max_val<<'\n';
    std::cout<<"Error (precond)  = "<<err_precond/max_val<<'\n';
    std::cout<<V0<<V1;

    for (Long j = 0; j < Ndisc; j++) {
      X[j*COORD_DIM+0] += V1[j*3+0] * dt;
      X[j*COORD_DIM+1] += V1[j*3+1] * dt;
    }
  //}

  return 0;
}


