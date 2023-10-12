#include <bie2d.hpp>
using namespace sctl;

const char data_path[]="data2/";
constexpr Integer COORD_DIM = 2;
constexpr Integer ElemOrder = 16;


template <class Real, class Disc> void RigidVelocityBasis(Matrix<Real>& V, const Vector<Disc>& disc_lst) {
  const Long N = NodeCount(disc_lst);
  if (V.Dim(0) != 3 || V.Dim(1) != N*COORD_DIM) V.ReInit(3, N*COORD_DIM);

  Vector<Real> X;
  Long offset = 0;
  for (Long k = 0; k < disc_lst.Dim(); k++) {
    const Real xx = disc_lst[k].Coord(0);
    const Real yy = disc_lst[k].Coord(1);
    const Real R = disc_lst[k].Radius();

    disc_lst[k].GetGeom(&X);
    const Long N_ = X.Dim()/COORD_DIM;

    const Real Rinv = 1/R;
    for (Long i = 0; i < N_; i++) {
      V[0][(offset+i)*COORD_DIM+0] = 1;
      V[0][(offset+i)*COORD_DIM+1] = 0;

      V[1][(offset+i)*COORD_DIM+0] = 0;
      V[1][(offset+i)*COORD_DIM+1] = 1;

      V[2][(offset+i)*COORD_DIM+0] =-(X[i*COORD_DIM+1]-yy) * Rinv;
      V[2][(offset+i)*COORD_DIM+1] = (X[i*COORD_DIM+0]-xx) * Rinv;
    }
    offset += N_;
  }
}

template <class Real, class PanelLstType> void RigidVelocityBasis(Matrix<Real>& V, const PanelLstType& panel_lst) {
  const Long N = panel_lst.Size() * panel_lst.ElemOrder();
  if (V.Dim(0) != 3 || V.Dim(1) != N*COORD_DIM) V.ReInit(3, N*COORD_DIM);

  Vector<Real> X;
  Long offset = 0;
  for (Long k = 0; k < panel_lst.DiscCount(); k++) {
    const Real xx = std::get<0>(panel_lst.DiscCoord(k));
    const Real yy = std::get<1>(panel_lst.DiscCoord(k));
    const Real R = panel_lst.DiscRadius();

    const auto& X = panel_lst.SurfCoord(k);
    const Long N_ = X.Dim() / COORD_DIM;

    const Real Rinv = 1/R;
    for (Long i = 0; i < N_; i++) {
      V[0][(offset+i)*COORD_DIM+0] = 1;
      V[0][(offset+i)*COORD_DIM+1] = 0;

      V[1][(offset+i)*COORD_DIM+0] = 0;
      V[1][(offset+i)*COORD_DIM+1] = 1;

      V[2][(offset+i)*COORD_DIM+0] =-(X[i*COORD_DIM+1]-yy) * Rinv;
      V[2][(offset+i)*COORD_DIM+1] = (X[i*COORD_DIM+0]-xx) * Rinv;
    }
    offset += N_;
  }
}


template <class Real, class Disc> void StokesOpMatrix(Matrix<Real>& K, Matrix<Real>& Msl, const Vector<Disc>& disc_lst, const Real tol) {
  using KerSL = Stokes2D_FxU;
  using KerDL = Stokes2D_DxU;

  const Long Ndisc = disc_lst.Dim();

  Vector<Real> Xt;
  GetGeom<Real>(disc_lst, &Xt);

  Matrix<Real> Mdl;
  LayerPotentialMatrix<KerSL>(Msl, disc_lst, Xt, tol);
  LayerPotentialMatrix<KerDL>(Mdl, disc_lst, Xt, tol);

  K = Mdl;
  for (Long i = 0; i < K.Dim(0); i++) K[i][i] += 0.5;

  { // add projection to ridig motion velocity
    Matrix<Real> V0;
    RigidVelocityBasis(V0, disc_lst);
    SCTL_ASSERT(V0.Dim(0)==3 && V0.Dim(1)==NodeCount(disc_lst)*COORD_DIM);

    Long offset = 0;
    Vector<Real> wts;
    for (Long k = 0; k < Ndisc; k++) {
      const Real inv_circumf = 1/(2*const_pi<Real>()*disc_lst[k].Radius());
      disc_lst[k].GetGeom(nullptr, nullptr, &wts);
      const Long N_ = wts.Dim() * COORD_DIM;
      for (Long i = 0; i < N_; i++) {
        const Long ii = (i>>1);
        for (Long j = 0; j < N_; j++) {
          for (Integer k = 0; k < 3; k++) {
            K[offset+i][offset+j] += inv_circumf * wts[ii] * V0[k][offset+i] * V0[k][offset+j];
          }
        }
      }
      offset += N_;
    }
  }
}

template <class Real> Vector<Real> DiscMobilitySolve_(const Vector<Real>& F, const Vector<Real>& X, const Real R, const Real tol = 1e-10) {
  const Long Nunif = 32; ////////////
  using DiscType = Disc<Real,ElemOrder>;
  const Long Ndisc = X.Dim() / COORD_DIM;

  //const auto disc_lst0 = GetAdapDisc<Real>(R, eps, 3);
  Vector<DiscType> disc_lst(Ndisc);
  for (Long i = 0; i < Ndisc; i++) {
    disc_lst[i] = DiscType(X[i*2+0], X[i*2+1], R, Nunif);
  }
  const Long N = NodeCount(disc_lst)*COORD_DIM;

  Matrix<Real> K, S; // Stokes mobility and single-layer operators
  StokesOpMatrix(K, S, disc_lst, tol);

  Matrix<Real> V0;
  RigidVelocityBasis(V0, disc_lst);
  SCTL_ASSERT(V0.Dim(0)==3 && V0.Dim(1)==N);

  Matrix<Real> nu; // boundary force
  { // Set nu
    nu.ReInit(1,N);
    Long offset = 0;
    for (Long i = 0; i < Ndisc; i++) {
      const Long N_ = disc_lst[i].NodeCount()*COORD_DIM;
      const Real inv_circumf = 1/(2*const_pi<Real>()*disc_lst[i].Radius());
      for (Long j = 0; j < N_; j++) {
        nu[0][offset+j] = 0;
        for (Integer k = 0; k < 3; k++) {
          nu[0][offset+j] += F[i*3+k] * V0[k][offset+j] * inv_circumf;
        }
      }
      offset += N_;
    }
  }

  const auto U0 = nu * S;
  const auto sigma = U0 * Matrix<Real>(K).pinv(); // Solve mobility

  Vector<Real> V;
  { // recover translation and rotation
    const Long Ndisc = disc_lst.Dim();
    V.ReInit(Ndisc*3);
    V = 0;

    Matrix<Real> V0;
    RigidVelocityBasis(V0, disc_lst);

    Long offset = 0;
    Vector<Real> wts;
    for (Long j = 0; j < Ndisc; j++) {
      disc_lst[j].GetGeom(nullptr, nullptr, &wts);
      const Long N_ = wts.Dim() * COORD_DIM;
      const Real inv_circumf = 1/(2*const_pi<Real>()*disc_lst[j].Radius());
      for (Long i = 0; i < N_; i++) {
        const Long ii = (i>>1);
        for (Integer k = 0; k < 3; k++) {
          V[j*3+k] += inv_circumf * wts[ii] * V0[k][offset+i] * sigma[0][offset+i];
        }
      }
      offset += N_;
    }
  }
  return V;
}

template <class Real> Vector<Real> DiscMobilitySolve(const Vector<Real>& F, const Vector<Real>& X, const Real R, const Real tol = 1e-10) {
  using PanelLstType = DiscPanelLst<Real,ElemOrder>;
  const Comm comm = Comm::Self();

  PanelLstType panel_lst(comm);
  panel_lst.Init(X, R, true);
  const Long N = panel_lst.Size() * ElemOrder * COORD_DIM;

  Stokes2D_FxU StokesSL_Ker;
  Stokes2D_DxU StokesDL_Ker;
  BoundaryIntegralOp<Real,Stokes2D_FxU> StokesSL_BIOp(StokesSL_Ker, false, comm);
  BoundaryIntegralOp<Real,Stokes2D_DxU> StokesDL_BIOp(StokesDL_Ker, false, comm);
  StokesSL_BIOp.AddElemList(panel_lst);
  StokesDL_BIOp.AddElemList(panel_lst);

  Matrix<Real> V0;
  RigidVelocityBasis(V0, panel_lst);
  SCTL_ASSERT(V0.Dim(0)==3 && V0.Dim(1)==N);

  const auto StokesMobilOp = [&StokesDL_BIOp,&panel_lst,&V0,&R](Vector<Real>* U, const Vector<Real> sigma) {
    StokesDL_BIOp.ComputePotential(*U, sigma);
    (*U) += 0.5 * sigma;

    Long offset = 0;
    for (Long k = 0; k < panel_lst.DiscCount(); k++) {
      const Real inv_circumf = 1/(2*const_pi<Real>()*R);
      const auto wts = panel_lst.SurfWts(k);
      const Long N = wts.Dim();

      Vector<Real> I(3); I = 0;
      for (Long i = 0; i < N; i++) {
        for (Long j = 0; j < 3; j++) {
          I[j] += wts[i] * sigma[offset+i*COORD_DIM+0] * V0[j][offset+i*COORD_DIM+0];
          I[j] += wts[i] * sigma[offset+i*COORD_DIM+1] * V0[j][offset+i*COORD_DIM+1];
        }
      }
      I *= inv_circumf;
      for (Long i = 0; i < N*COORD_DIM; i++) {
        for (Long j = 0; j < 3; j++) {
          (*U)[offset+i] += V0[j][offset+i] * I[j];
        }
      }
      offset += N*COORD_DIM;
    }
  };

  Vector<Real> nu; // boundary force
  { // Set nu
    nu.ReInit(N);
    Long offset = 0;
    for (Long i = 0; i < panel_lst.DiscCount(); i++) {
      const Long N_ = panel_lst.SurfCoord(i).Dim();
      const Real inv_circumf = 1/(2*const_pi<Real>()*R);
      for (Long j = 0; j < N_; j++) {
        nu[offset+j] = 0;
        for (Integer k = 0; k < 3; k++) {
          nu[offset+j] += F[i*3+k] * V0[k][offset+j] * inv_circumf;
        }
      }
      offset += N_;
    }
  }

  Vector<Real> U0;
  StokesSL_BIOp.ComputePotential(U0, nu);

  Vector<Real> sigma;
  ParallelSolver<Real> solver;
  solver(&sigma, StokesMobilOp, U0, 1e-8, 100);

  Vector<Real> V;
  { // recover translation and rotation
    const Long Ndisc = panel_lst.DiscCount();
    V.ReInit(Ndisc*3);
    V = 0;

    //Matrix<Real> V0;
    //RigidVelocityBasis(V0, disc_lst);

    Long offset = 0;
    Vector<Real> wts;
    for (Long j = 0; j < Ndisc; j++) {
      const auto wts = panel_lst.SurfWts(j);
      const Long N_ = wts.Dim() * COORD_DIM;
      const Real inv_circumf = 1/(2*const_pi<Real>()*R);
      for (Long i = 0; i < N_; i++) {
        const Long ii = (i>>1);
        for (Integer k = 0; k < 3; k++) {
          V[j*3+k] += inv_circumf * wts[ii] * V0[k][offset+i] * sigma[offset+i];
        }
      }
      offset += N_;
    }
  }
  return V;
}

int main() {
  using Real = double;
  const Real tol = machine_eps<Real>();

  const Real dt = 0.1;
  const Real R = 0.75;
  const Long Ndisc = 2;

  Vector<Real> X(Ndisc*COORD_DIM), F(Ndisc*3);
  { // Set F, X
    X = 0;
    X[0] = -0.8; X[1] = 0.0;
    X[2] =  0.8; X[3] = 0.0;
    //X[4] =  0.0; X[5] =-0.5;

    F = 0;
    F[0] = 0; F[1] = 1; F[2] = 0;
    F[3] = 0; F[4] = 0; F[5] = 0;
    //F[6] = 0; F[7] = 0; F[8] = 0;
  }

  //for (Long i = 0; i < 1000; i++) {
    const Vector<Real> V = DiscMobilitySolve(F, X, R, tol);
    std::cout<<X;
    std::cout<<V;
    for (Long j = 0; j < Ndisc; j++) {
      X[j*COORD_DIM+0] += V[j*3+0] * dt;
      X[j*COORD_DIM+1] += V[j*3+1] * dt;
    }
  //}

  return 0;
}


