#include <bie2d.hpp>
using namespace sctl;

constexpr Integer ElemOrder = 16;
using RefReal = long double; // reference solution precision
using Real = double;

const char data_path[]="data2/";
constexpr Integer COORD_DIM = 2;

template <Integer ElemOrder, class Real> Vector<Real> DiscMobilitySolve(const Vector<Real>& F, const Vector<Real>& X, const Real R, const ICIPType icip_type, const Real tol) {
  const Comm comm = Comm::Self();
  DiscMobility<Real,ElemOrder> disc_mobility(comm);

  Vector<Real> V;
  disc_mobility.Init(X, R, tol, icip_type);
  disc_mobility.Solve(V, F, Vector<Real>(), 0.0000001, 2000);
  return V;
}

int main() {
  const Real tol = machine_eps<Real>()*64;

  const Real dt = 0.1;
  const Real R = 0.75;
  const Long Ndisc = 2;

  Vector<Real> X(Ndisc*COORD_DIM), F(Ndisc*3);
  { // Set F, X
    const Real eps = 4.6509e-9;

    X = 0;
    // X[0] = -(R+eps/2)*sqrt<Real>(0.5); X[1] = -(R+eps/2)*sqrt<Real>(0.5);
    // X[2] =  (R+eps/2)*sqrt<Real>(0.5); X[3] =  (R+eps/2)*sqrt<Real>(0.5);
    //X[4] =  0.0; X[5] =-0.5;
    X[0] = -R/2*eps - R; X[1] = 0;
    X[2] =  R/2*eps + R; X[3] =  0;

    F = 0;
    F[0] = 0; F[1] = 1; F[2] = 0;
    F[3] = 0; F[4] = 0; F[5] = 0;
    //F[6] = 0; F[7] = 0; F[8] = 0;
  }

  const auto real2quad = [](const Vector<Real>& v) {
    Vector<RefReal> w(v.Dim());
    for (Long i = 0; i < v.Dim(); i++) w[i] = (RefReal)v[i];
    return w;
  };
  const auto quad2real = [](const Vector<RefReal>& v) {
    Vector<Real> w(v.Dim());
    for (Long i = 0; i < v.Dim(); i++) w[i] = (Real)v[i];
    return w;
  };

  //for (Long i = 0; i < 1000; i++) {
    const auto V0 = quad2real(DiscMobilitySolve<ElemOrder+8,RefReal>(real2quad(F), real2quad(X), R, ICIPType::Adaptive, tol*1e-3));
    const auto V1 = DiscMobilitySolve<ElemOrder>(F, X, R, ICIPType::Adaptive, tol);
    const auto V2 = DiscMobilitySolve<ElemOrder>(F, X, R, ICIPType::Compress, tol);
    const auto V3 = DiscMobilitySolve<ElemOrder>(F, X, R, ICIPType::Precond , tol);

    Real err_adap = 0, err_comp = 0, err_prec = 0, max_val = 0;
    for (const auto x : V0) max_val = std::max<Real>(max_val, fabs(x));
    for (const auto x : V0-V1) err_adap = std::max<Real>(err_adap, fabs(x));
    for (const auto x : V0-V2) err_comp = std::max<Real>(err_comp, fabs(x));
    for (const auto x : V0-V3) err_prec = std::max<Real>(err_prec, fabs(x));
    std::cout<<"Error (adaptive) = "<<err_adap/max_val<<'\n';
    std::cout<<"Error (compress) = "<<err_comp/max_val<<'\n';
    std::cout<<"Error (precond)  = "<<err_prec/max_val<<'\n';
    std::cout<<V0;

    for (Long j = 0; j < Ndisc; j++) {
      X[j*COORD_DIM+0] += V1[j*3+0] * dt;
      X[j*COORD_DIM+1] += V1[j*3+1] * dt;
    }
  //}

  return 0;
}


