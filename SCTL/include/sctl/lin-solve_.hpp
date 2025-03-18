#ifndef _SCTL_LIN_SOLVE_HPP_
#define _SCTL_LIN_SOLVE_HPP_

#include <sctl/common.hpp>
#include SCTL_INCLUDE(comm.hpp)
#include SCTL_INCLUDE(mem_mgr.hpp)
#include SCTL_INCLUDE(math_utils.hpp)

#include <functional>
#include <list>

namespace SCTL_NAMESPACE {

template <class ValueType> class Vector;
template <class ValueType> class Matrix;

#if 1
template <class Real> class KrylovPrecond {
  public:

    KrylovPrecond() : N_(0) {}

    Long Size() const { return N_; }

    Long Rank() const {
      Long rank = 0;
      for (auto it = mat_lst.begin(); it != mat_lst.end(); it++) {
        rank += it->Dim(1);
        it++;
      }
      return rank;
    }

    void RankOneUpdate(const Vector<Real>& Ax, const Vector<Real>& x) {} /////////////////////////
    void SetRHS(const Vector<Real>& b) {} //////////////////////////
    template <class T> void AugmentSubspace(const T& p) {} //////////////////////////

    void Append(const Matrix<Real>& Qt, const Matrix<Real>& U) {
      SCTL_ASSERT(Qt.Dim(0) == U.Dim(1));
      SCTL_ASSERT(Qt.Dim(1) == U.Dim(0));
      if (Qt.Dim(0) != N_) { // clear
        mat_lst.clear();
        N_ = Qt.Dim(0);
      }

      mat_lst.push_front(U);
      mat_lst.push_front(Qt);

      if (mat_lst.size() > 57*2+2) { ///////////////////// prune
        const Long N = mat_lst.size() / 2;

        Long iter_cnt = 0;
        std::list<Matrix<Real>> mat_lst_trunc;
        for (auto it = mat_lst.begin(); it != mat_lst.end(); it++) {
          const auto& Qt = *it;
          it++;
          const auto& U = *it;

          if (iter_cnt < 1 || iter_cnt > N-8) {
            mat_lst_trunc.push_back(Qt);
            mat_lst_trunc.push_back(U);
          }
          iter_cnt++;
        }
        mat_lst.swap(mat_lst_trunc);
      }
    }

    void Apply(Vector<Real>& y) const {
      if (N_ != y.Dim()) return;

      Matrix<Real> y_(1, N_, y.begin(), false);
      for (auto it = mat_lst.begin(); it != mat_lst.end(); it++) {
        const auto& Qt = *it;
        it++;
        const auto& U = *it;
        y_ += (y_ * Qt) * U;
      }
    }

  private:

    Long N_;
    std::list<Matrix<Real>> mat_lst;
};
#endif

#if 0
template <class Real> class KrylovPrecond {
  public:

    KrylovPrecond(const Real tol = 1e-16, const Long length = 57) : N_(-1), tol2(tol*tol), length_(length) {}

    Long Size() const { return N_; }

    Long Rank() const { return Qt.Dim(1); }

    void RankOneUpdate(const Vector<Real>& Ax, const Vector<Real>& x) {
      if (N_ < 0) N_ = x.Dim();
      SCTL_ASSERT(x.Dim() == N_);
      SCTL_ASSERT(Ax.Dim() == N_);

      const auto norm2 = [](const Vector<Real>& x){
        Real norm2 = 0;
        for (const auto& a : x) norm2 += a*a;
        return sqrt<Real>(norm2);
      };
      const Real norm2_x = norm2(x);
      const Real norm2_Ax = norm2(Ax);
      if (norm2_x == 0 || norm2_Ax == 0) return;

      const Real scal = 1/norm2_x;
      for (const auto& a : x) x_.PushBack(a*scal);
      for (const auto& a : Ax) Ax_.PushBack(a*scal);

      if (Ax_vec.Dim() == 0) {
        Ax_vec.ReInit(1);
        x_vec.ReInit(1);
      }
      Ax_vec[0].PushBack(Ax); /////////////////
      x_vec[0].PushBack(x); /////////////////
    }
    void SetRHS(const Vector<Real>& b) {
      rhs.PushBack(b);
    }

    void AugmentSubspace(const KrylovPrecond& p) {
      if (p.Size() != Size()) {
        //static Long rank = 0;
        //if (Rank() > rank) {
        //  rank = Rank();
        //  Qt.Write("Qt.mat");
        //  W.Write("W.mat");
        //}
        (*this) = p;
        return;
      }
      if (Ax_vec.Dim() > length_) {
        KrylovPrecond<Real> pp;
        for (Long i = Ax_vec.Dim()-1; i >= 0; i--) {
          if (0 <= i && i < Ax_vec.Dim() - 7) continue;
          KrylovPrecond<Real> p(1e-16, 1000);
          for (Long j = 0; j < Ax_vec[i].Dim(); j++) {
            p.RankOneUpdate(Ax_vec[i][j], x_vec[i][j]);
          }
          if (pp.Size() != p.Size()) pp = p;
          else pp.AugmentSubspace(p);
        }
        //if (pp.Size() != p.Size()) pp = p;
        //else pp.AugmentSubspace(p);
        (*this) = pp;
      }

      // (I-Qt2*Q2) * (I + Qt1*W1) + Qt2*(W2+Q2)
      // => I + Qt1*W1 + Qt2*(W2 - Q2*Qt1*W1)

      p.Setup();
      Setup();

      auto Qt1 = Qt;
      auto W1 = W;

      auto Qt2 = p.Qt;
      auto W2 = p.W - (Qt2.Transpose()*Qt1) * W1;

      //// orthogonalize Q1 wrt Q2
      //const auto Q2_Qt1 = Qt2.Transpose() * Qt1;
      //Qt1 -= Qt2 * Q2_Qt1;
      //W2 += Q2_Qt1 * W1;

      //{ // Orthogonalize Q1
      //  const auto Q1 = Qt1.Transpose();
      //  const auto Ainv_Q1 = Q1 + (Q1 * Qt1) * W1;

      //  KrylovPrecond<Real> pp;
      //  Vector<Real> Ax(N_), x(N_);
      //  for (Long k = 0; k < Qt1.Dim(1); k++) {
      //    for (Long i = 0; i < N_; i++) {
      //      Ax[i] = Qt1[i][k];
      //      x[i]  = Ainv_Q1[k][i];
      //    }
      //    pp.RankOneUpdate(Ax, x);
      //  }
      //  pp.Setup();
      //  Qt1 = pp.Qt;
      //  W1 = pp.W;
      //}

      if (1) { /////////////////////////////////////////////////////////////
        Vector<Vector<Vector<Real>>> Ax_vec_, x_vec_;
        for (const auto& a : p.Ax_vec) Ax_vec_.PushBack(a);
        for (const auto& a :   Ax_vec) Ax_vec_.PushBack(a);
        Ax_vec = Ax_vec_;

        for (const auto& a : p.x_vec) x_vec_.PushBack(a);
        for (const auto& a :   x_vec) x_vec_.PushBack(a);
        x_vec = x_vec_;

        for (const auto& v : p.rhs) rhs.PushBack(v);
      }
      if (0) { ///////////////////////////
        Vector<Matrix<Real>> Q_vec(Ax_vec.Dim());
        for (Long i = 0 ; i < Ax_vec.Dim(); i++) {
          const Long L = Ax_vec[i].Dim();
          Q_vec[i].ReInit(L, N_);
          for (Long j = 0; j < L; j++) { // Set Q_vec[i]
            for (Long k = 0; k < N_; k++) {
              Q_vec[i][j][k] = Ax_vec[i][j][k];
            }
          }
          for (Long j = 0; j < i; j++) { // orthogonalize Q_vec[i] with Q_vec[0:i-1]
            Q_vec[i] -= (Q_vec[i] * Q_vec[j].Transpose()) * Q_vec[j];
          }
          for (Long j0 = 0; j0 < L; j0++) {
            for (Long j1 = 0; j1 < j0; j1++) { // orthogonalize
              Real inner_prod = 0;
              for (Long k = 0; k < N_; k++) inner_prod += Q_vec[i][j0][k] * Q_vec[i][j1][k];
              for (Long k = 0; k < N_; k++) Q_vec[i][j0][k] -= inner_prod * Q_vec[i][j1][k];
            }
            { // normalize
              Real inner_prod = 0;
              for (Long k = 0; k < N_; k++) inner_prod += Q_vec[i][j0][k] * Q_vec[i][j0][k];
              for (Long k = 0; k < N_; k++) Q_vec[i][j0][k] *= 1/sqrt<Real>(inner_prod);
            }
          }
        }

        for (Long i = 0; i < rhs.Dim(); i++) {
          Matrix<Real> b(1, N_, rhs[i].begin());
          for (Long j = 0; j < Ax_vec.Dim(); j++) {
            b -= (b * Q_vec[j].Transpose()) * Q_vec[j];
            printf(" %.2e ", sqrt<Real>((b*b.Transpose())[0][0]));
          }
          std::cout<<'\n';
        }
      }



      const Long M1 = Qt1.Dim(1);
      const Long M2 = Qt2.Dim(1);
      W.ReInit(M1+M2, N_);
      Matrix<Real> Q(M1+M2, N_);
      for (Long j = 0; j < M1; j++) {
        for (Long i = 0; i < N_; i++) Q[j][i] = Qt1[i][j];
        for (Long i = 0; i < N_; i++) W[j][i] = W1[j][i];
      }
      for (Long j = 0; j < M2; j++) {
        for (Long i = 0; i < N_; i++) Q[M1+j][i] = Qt2[i][j];
        for (Long i = 0; i < N_; i++) W[M1+j][i] = W2[j][i];
      }


      //for (Long k = 0; k < M1+M2; k++) {
      //  Real norm2_Q = 0, norm2_W = 0;
      //  for (Long i = 0; i < N_; i++) norm2_Q += Q[k][i] * Q[k][i];
      //  for (Long i = 0; i < N_; i++) norm2_W += (W[k][i]+Q[k][i]) * (W[k][i]+Q[k][i]);
      //  if (norm2_Q * norm2_W < tol2) {
      //    for (Long i = 0; i < N_; i++) Q[k][i] = 0;
      //    for (Long i = 0; i < N_; i++) W[k][i] = 0;
      //  }
      //  if (norm2_Q * norm2_W > 1/tol2) {
      //    for (Long i = 0; i < N_; i++) Q[k][i] = 0;
      //    for (Long i = 0; i < N_; i++) W[k][i] = 0;
      //  }
      //}

      Qt = Q.Transpose();
      //Compress();


      if (0) { ///////////////////////////////////////////////////////////////////
        auto E = Q * Qt;
        for (Long i = 0; i < E.Dim(0); i++) if (E[i][i] != 0) E[i][i] -= 1;
        Real max_err = 0;
        for (const auto a : E) max_err = std::max<Real>(max_err, fabs(a));

        Long rank = 0;
        Matrix<Real> U,S,Vt, W_ = W+Q;
        W_.SVD(U,S,Vt);
        Real max_val = S[0][0], min_val = S[0][0];
        for (Long i = 0; i < std::min(S.Dim(0), S.Dim(1)); i++) {
          max_val = std::max<Real>(max_val, fabs(S[i][i]));
          if (S[i][i] != 0) {
            rank++;
            min_val = std::min<Real>(min_val, fabs(S[i][i]));
          }
        }

        std::cout<<"GMRES rank = "<<rank<<"    Q-error = "<<max_err<<"    Smax/Smin = "<<max_val<<" / "<<min_val<<'\n';
      }


      // TODO: accumulate all previous Ax, x and see how whell they agree with the preconditioner
      // Or see if we still run into issue if we construct the preconditioner from all previous Ax, x
      if (0) if (Ax_vec.Dim() > length_) {
        Vector<Real> x;
        KrylovPrecond<Real> pp;
        for (Long i = Ax_vec.Dim()-1; i >= 0; i--) {
          if (0 <= i && i <= Ax_vec.Dim()-1 - 7) continue;
          for (Long j = 0; j < Ax_vec[i].Dim(); j++) {
            x = Ax_vec[i][j];
            Apply(x);
            pp.RankOneUpdate(Ax_vec[i][j], x);
          }
        }
        (*this) = pp;
      }
    }

    void AugmentSubspace1(const KrylovPrecond& p) {
      if (p.Size() != Size()) {
        (*this) = p;
        return;
      }
      p.Setup();
      Setup();

      auto Qt1 = Qt;
      auto W1 = W;

      auto Qt2 = p.Qt;
      auto W2 = p.W + (p.W * Qt) * W;


      const Long M1 = Qt1.Dim(1);
      const Long M2 = Qt2.Dim(1);
      W.ReInit(M1+M2, N_);
      Matrix<Real> Q(M1+M2, N_);
      for (Long i = 0; i < N_; i++) {
        for (Long j = 0; j < M1; j++) Q[   j][i] = Qt1[i][j];
        for (Long j = 0; j < M2; j++) Q[M1+j][i] = Qt2[i][j];
      }
      for (Long j = 0; j < M1; j++) {
        for (Long i = 0; i < N_; i++) W[   j][i] = W1[j][i];
      }
      for (Long j = 0; j < M2; j++) {
        for (Long i = 0; i < N_; i++) W[M1+j][i] = W2[j][i];
      }



      //for (Long j = 0; j < i; j++) {
      //  for (Long i = M1; i < M1+M2; i++) {
      //    Real inner_prod = 0;
      //    for (Long k = 0; k < N_; k++) inner_prod += Q[i][k] * Q[j][k];
      //    for (Long k = 0; k < N_; k++) Q[i][k] -= inner_prod * Q[j][k];
      //    for (Long k = 0; k < N_; k++) W[j][k] += inner_prod * W[i][k];

      //    inner_prod = 0;
      //    for (Long k = 0; k < N_; k++) inner_prod += Q[i][k] * Q[i][k];
      //    const Real norm_qi = sqrt<Real>(inner_prod);
      //    const Real inv_norm_qi = 1 / norm_qi;
      //    for (Long k = 0; k < N_; k++) Q[i][k] *= inv_norm_qi;
      //    for (Long k = 0; k < N_; k++) W[i][k] *= norm_qi;
      //  }
      //}

      for (Long i = M1; i < M1+M2; i++) {
        for (Long j = 0; j < i; j++) {
          Real inner_prod = 0;
          for (Long k = 0; k < N_; k++) inner_prod += Q[i][k] * Q[j][k];
          for (Long k = 0; k < N_; k++) Q[i][k] -= inner_prod * Q[j][k];
          for (Long k = 0; k < N_; k++) W[j][k] += inner_prod * W[i][k];
        }
        { // normalize Q[i]
          Real inner_prod = 0;
          for (Long k = 0; k < N_; k++) inner_prod += Q[i][k] * Q[i][k];
          const Real norm_qi = sqrt<Real>(inner_prod);
          const Real inv_norm_qi = 1 / norm_qi;
          for (Long k = 0; k < N_; k++) Q[i][k] *= inv_norm_qi;
          for (Long k = 0; k < N_; k++) W[i][k] *= norm_qi;

          inner_prod = 0;
          for (Long k = 0; k < N_; k++) inner_prod += W[i][k] * W[i][k];
          if (inner_prod < tol2) {
            for (Long k = 0; k < N_; k++) Q[i][k] = 0;
            for (Long k = 0; k < N_; k++) W[i][k] = 0;
          }
        }
      }

      Qt = Q.Transpose();
      //if (Rank() > 2500) Compress();
    }

    void Append(const Matrix<Real>& Qt, const Matrix<Real>& U) {} ///////////////////////////////

    void Apply(Vector<Real>& y) const {
      if (N_ != y.Dim()) return;
      Setup();

      Matrix<Real> y_(1, N_, y.begin(), false);
      y_ = y_ + (y_ * Qt) * W;
    }

  private:

    void Setup() const {
      if (x_.Dim() == 0) return;
      const Long K = x_.Dim() / N_;
      Matrix<Real> Q_(K, N_); Q_ = 0;
      Matrix<Real> W_(K, N_); W_ = 0;

      auto x = x_;
      auto Ax = Ax_;
      Vector<Real> norm2_Ax(K);
      for (Long j = 0; j < K; j++) {
        norm2_Ax[j] = 0;
        for (Long i = 0; i < N_; i++) {
          norm2_Ax[j] += Ax[j*N_+i] * Ax[j*N_+i];
        }
      }
      //std::cout<<"norm2_Ax = "<<norm2_Ax<<'\n';

      Long rank = K;
      for (Long k = 0; k < K; k++) {
        Real scal = 1/sqrt<Real>(norm2_Ax[k]);
        for (Long i = 0; i < N_; i++) { // Q_[k] <-- Ax[k] * scal
          Q_[k][i] = Ax[k*N_+i] * scal;
          W_[k][i] = x[k*N_+i] * scal;
        }

        #pragma omp parallel for schedule(static)
        for (Long j = 0; j < K; j++) { // orthogonalize Ax[k+1:K] with Q_[k];
          if (norm2_Ax[j] > 0) {
            Real inner_prod  = 0,  norm2_Ax_ = 0;
            for (Long i = 0; i < N_; i++) {
              inner_prod += Q_[k][i] * Ax[j*N_+i];
            }
            for (Long i = 0; i < N_; i++) {
              Ax[j*N_+i] -= Q_[k][i] * inner_prod;
              x[j*N_+i] -= W_[k][i] * inner_prod;
              norm2_Ax_ += Ax[j*N_+i] * Ax[j*N_+i];
            }
            norm2_Ax[j] = norm2_Ax_;
          }
        }
      }

      { // Q <-- Q_[1:rank], W <-- W_[1:rank] - Q
        W.ReInit(rank, N_);
        Qt.ReInit(N_, rank);
        for (Long i = 0; i < rank; i++) {
          for (Long j = 0; j < N_; j++) {
            W[i][j] = W_[i][j] - Q_[i][j];
            Qt[j][i] = Q_[i][j];
          }
        }
      }

      Ax_.ReInit(0);
      x_.ReInit(0);
    }

    void Compress() {
      // TODO: do it using QR of Q matrix
      const auto Q = Qt.Transpose();
      const Matrix<Real> Ainv_Q = Q + (Q * Qt) * W;
      for (Long k = 0; k < Q.Dim(0); k++) {
        Vector<Real> Ax(N_, (Iterator<Real>)     Q.begin() + k*N_, false);
        Vector<Real>  x(N_, (Iterator<Real>)Ainv_Q.begin() + k*N_, false);
        RankOneUpdate(Ax, x);
      }
    }

    Long N_;
    Real tol2;
    Long length_;

    mutable Matrix<Real> Qt, W;
    mutable Vector<Real> Ax_, x_;
    //mutable Vector<Real> Ax__, x__;

    mutable Vector<Vector<Real>> rhs;
    mutable Vector<Vector<Vector<Real>>> Ax_vec, x_vec;
};
#endif

#if 0
template <class Real> class KrylovPrecond {
  public:

    KrylovPrecond(const Real tol = 1e-14) : setup_status(false), N_(-1), tol_(tol) {}

    Long Size() const { return N_; }

    Long Rank() const { return Q.Dim(0); }

    const Vector<Real>& Get_Ax() { return Ax_; }
    const Vector<Real>& Get_x() { return x_; }

    void RankOneUpdate(const Vector<Real>& Ax, const Vector<Real>& x) {
      setup_status = false;
      if (N_ < 0) N_ = x.Dim();
      SCTL_ASSERT(x.Dim() == N_);
      SCTL_ASSERT(Ax.Dim() == N_);
      const Real scal = [&Ax](){
        Real norm2 = 0;
        for (const auto& a : Ax) norm2 += a*a;
        std::cout<<sqrt<Real>(norm2)<<'\n';
        return 1/sqrt<Real>(norm2);
      }();
      for (const auto& a : x) x_.PushBack(a*scal);
      for (const auto& a : Ax) Ax_.PushBack(a*scal);
    }

    void AugmentSubspace(const KrylovPrecond& p) {
      if (p.Size() != Size()) return;

      this->Setup();
      Ax_.ReInit(Q.Dim(0)*Q.Dim(1), Q.begin());
      x_.ReInit(Ainv_Q.Dim(0)*Ainv_Q.Dim(1), Ainv_Q.begin());

      Matrix<Real> Ax0(p.Ax_.Dim()/N_, N_, (Iterator<Real>)p.Ax_.begin());
      Matrix<Real> x0(p.x_.Dim()/N_, N_, (Iterator<Real>)p.x_.begin());

      const auto Ax0_Qt = Ax0 * Q.Transpose();
      Ax0 -= Ax0_Qt * Q;
      x0 -= Ax0_Qt * Ainv_Q;

      for (const auto& a : Ax0) Ax_.PushBack(a);
      for (const auto& a : x0) x_.PushBack(a);
      setup_status = false;
    }

    void ApplyInverse(Vector<Real>& y) const {
      if (N_ != y.Dim()) return;
      Setup();

      Matrix<Real> y_(1, N_, y.begin(), false);
      const auto y_Qt = (y_ * Q.Transpose());

      y_ = (y_ - y_Qt * Q) + y_Qt * Ainv_Q;
    }

  private:

    void Setup() const {
      if (setup_status) return;
      const Long K = x_.Dim() / N_;
      Ainv_Q.ReInit(K, N_); Ainv_Q = 0;
      Q.ReInit(K, N_); Q = 0;

      auto x = x_;
      auto Ax = Ax_;
      Vector<Real> norm2_Ax(K);
      for (Long j = 0; j < K; j++) {
        norm2_Ax[j] = 0;
        for (Long i = 0; i < N_; i++) {
          norm2_Ax[j] += Ax[j*N_+i] * Ax[j*N_+i];
        }
      }

      Long rank = K;
      for (Long k = 0; k < K; k++) {
        Long pivot = 0;
        Real pivot_norm2 = norm2_Ax[0];
        for (Long j = 1; j < K; j++) { // find pivot
          if (norm2_Ax[j] > pivot_norm2) {
            pivot = j;
            pivot_norm2 = norm2_Ax[j];
          }
        }

        Real scal = (pivot_norm2 > tol_ ? 1/sqrt<Real>(pivot_norm2) : 0);
        if (scal == 0) {
          rank = k;
          break;
        }
        for (Long i = 0; i < N_; i++) { // Q[k] <-- Ax[pivot] * scal
          Q[k][i] = Ax[pivot*N_+i] * scal;
          Ainv_Q[k][i] = x[pivot*N_+i] * scal;
        }

        #pragma omp parallel for schedule(static)
        for (Long j = 0; j < K; j++) { // orthogonalize Ax[k+1:K] with Q[k];
          if (norm2_Ax[j] > 0) {
            norm2_Ax[j] = 0;
            Real inner_prod  = 0;
            for (Long i = 0; i < N_; i++) {
              inner_prod += Q[k][i] * Ax[j*N_+i];
            }
            for (Long i = 0; i < N_; i++) {
              Ax[j*N_+i] -= Q[k][i] * inner_prod;
              x[j*N_+i] -= Ainv_Q[k][i] * inner_prod;
              norm2_Ax[j] += Ax[j*N_+i] * Ax[j*N_+i];
            }
          }
        }
      }

      { // resize Q, Ainv_Q to rank x N_
        Matrix<Real> Q_(rank, N_, Ax.begin(), false);
        Matrix<Real> Ainv_Q_(rank, N_, x.begin(), false);
        Ainv_Q_ = Ainv_Q;
        Q_ = Q;
        Ainv_Q.ReInit(rank, N_, Ainv_Q_.begin());
        Q.ReInit(rank, N_, Q_.begin());
      }

      setup_status = true;
    }

    mutable bool setup_status;
    mutable Matrix<Real> Q, Ainv_Q;

    Long N_;
    Real tol_;
    Vector<Real> Ax_, x_;
};
#endif

template <class Real> class GMRES {

 public:
  using ParallelOp = std::function<void(Vector<Real>*, const Vector<Real>&)>;

  GMRES(const Comm& comm = Comm::Self(), bool verbose = true) : comm_(comm), verbose_(verbose) {}

  void operator()(Vector<Real>* x, const ParallelOp& A, const Vector<Real>& b, const Real tol, const Integer max_iter = -1, const bool use_abs_tol = false, Long* solve_iter=nullptr, KrylovPrecond<Real>* precond=nullptr);

  static void test(Long N = 15) {
    srand48(0);
    Matrix<Real> A(N, N);
    Vector<Real> b(N), x;
    for (Long i = 0; i < N; i++) {
      b[i] = drand48();
      for (Long j = 0; j < N; j++) {
        A[i][j] = drand48();
      }
    }
    auto LinOp = [&A](Vector<Real>* Ax, const Vector<Real>& x) {
      const Long N = x.Dim();
      Ax->ReInit(N);
      Matrix<Real> Ax_(N, 1, Ax->begin(), false);
      Ax_ = A * Matrix<Real>(N, 1, (Iterator<Real>)x.begin(), false);
    };

    Long solve_iter;
    GMRES<Real> solver;
    solver(&x, LinOp, b, 1e-10, -1, false, &solve_iter);

    auto print_error = [N,&A,&b](const Vector<Real>& x) {
      Real max_err = 0;
      auto Merr = A*Matrix<Real>(N, 1, (Iterator<Real>)x.begin(), false) - Matrix<Real>(N, 1, b.begin(), false);
      for (const auto& a : Merr) max_err = std::max(max_err, fabs(a));
      std::cout<<"Maximum error = "<<max_err<<'\n';
    };
    print_error(x);
    std::cout<<"GMRES iterations = "<<solve_iter<<'\n';
  }

 private:
  void GenericGMRES(Vector<Real>* x, const ParallelOp& A, const Vector<Real>& b, const Real tol, Integer max_iter, const bool use_abs_tol, Long* solve_iter, KrylovPrecond<Real>* precond);

  Comm comm_;
  bool verbose_;
};

}  // end namespace

namespace SCTL_NAMESPACE {

template <class Real> static Real inner_prod(const Vector<Real>& x, const Vector<Real>& y, const Comm& comm) {
  Real x_dot_y = 0;
  Long N = x.Dim();
  SCTL_ASSERT(y.Dim() == N);
  for (Long i = 0; i < N; i++) x_dot_y += x[i] * y[i];

  Real x_dot_y_glb = 0;
  comm.Allreduce(Ptr2ConstItr<Real>(&x_dot_y, 1), Ptr2Itr<Real>(&x_dot_y_glb, 1), 1, Comm::CommOp::SUM);

  return x_dot_y_glb;
}

template <class Real> inline void GMRES<Real>::GenericGMRES(Vector<Real>* x, const ParallelOp& A, const Vector<Real>& b, Real tol, Integer max_iter, bool use_abs_tol, Long* solve_iter, KrylovPrecond<Real>* precond) {
  const Long N = b.Dim();
  KrylovPrecond<Real> precond_;
  if (max_iter < 0) { // set max_iter
    StaticArray<Long,2> NN{N,0};
    comm_.Allreduce(NN+0, NN+1, 1, Comm::CommOp::SUM);
    max_iter = NN[1];
  }
  static constexpr Real ARRAY_RESIZE_FACTOR = 1.618;

  Vector<Real> Q_mat, H_mat;
  auto ResizeVector = [](Vector<Real>& v, const Long N0) {
    if (v.Dim() < N0) {
      Vector<Real> v_(N0);
      for (Long i = 0; i < v.Dim(); i++) v_[i] = v[i];
      for (Long i = v.Dim(); i < N0; i++) v_[i] = 0;
      v.Swap(v_);
    }
  };
  auto Q_row = [N,&Q_mat,&ResizeVector](Long i) -> Iterator<Real> {
    const Long idx = i*N;
    if (Q_mat.Dim() <= idx+N) {
      ResizeVector(Q_mat, (Long)((idx+N)*ARRAY_RESIZE_FACTOR));
    }
    return Q_mat.begin() + idx;
  };
  auto Q = [&Q_row](Long i, Long j) -> Real& {
    return Q_row(i)[j];
  };
  auto H_row = [&H_mat,&ResizeVector](Long i) -> Iterator<Real> {
    const Long idx = i*(i+1)/2;
    if (H_mat.Dim() <= idx+i+1) ResizeVector(H_mat, (Long)((idx+i+1)*ARRAY_RESIZE_FACTOR));
    return H_mat.begin() + idx;
  };
  auto H = [&H_row](Long i, Long j) -> Real& {
    return H_row(i)[j];
  };

  auto apply_givens_rotation = [](Vector<Real>& h, Real& cs_k, Real& sn_k, const Vector<Real>& cs, const Vector<Real>& sn, const Long k) {
    // apply for ith row
    for (Long i = 0; i < k; i++) {
      Real temp = cs[i] * h[i] + sn[i] * h[i+1];
      h[i+1]   = -sn[i] * h[i] + cs[i] * h[i+1];
      h[i]     = temp;
    }

    // update the next sin cos values for rotation
    const Real t = sqrt<Real>(h[k]*h[k] + h[k+1]*h[k+1]);
    cs_k = h[k] / t;
    sn_k = h[k+1] / t;

    // eliminate H(i + 1, i)
    h[k] = cs_k * h[k] + sn_k * h[k+1];
    h[k+1] = 0.0;
  };
  auto arnoldi = [this,N,&Q_row,&Q,&precond,&precond_](Vector<Real>& h, Vector<Real>& q, const ParallelOp& A, const Long k) {
    Vector<Real> q_k(N, Q_row(k), precond);
    if (precond) precond->Apply(q_k);
    A(&q, q_k);
    if (precond) precond_.RankOneUpdate(q, q_k);

    for (Long i = 0; i < k+1; i++) { // Modified Gram-Schmidt, keeping the Hessenberg matrix
      h[i] = inner_prod(q, Vector<Real>(N, Q_row(i), false), comm_);
      for (Long j = 0; j < N; j++) {
        q[j] -= h[i] * Q(i,j);
      }
    }
    for (Long i = 0; i < k+1; i++) { // re-orthogonalize (more stable)
      const Real h_ = inner_prod(q, Vector<Real>(N, Q_row(i), false), comm_);
      for (Long j = 0; j < N; j++) {
        q[j] -= h_ * Q(i,j);
      }
      h[i] += h_;
    }
    h[k+1] = sqrt<Real>(inner_prod(q, q, comm_));
    q *= 1/h[k+1];
  };

  Vector<Real> r;
  if (x->Dim() == N) { // r = b - A * x;
    Vector<Real> Ax;
    A(&Ax, *x);
    if (precond) precond_.RankOneUpdate(Ax, *x);
    r = b - Ax;
  } else {
    r = b;
    x->ReInit(N);
    x->SetZero();
  }

  const Real b_norm = sqrt<Real>(inner_prod(b, b, comm_));
  const Real abs_tol = tol * (use_abs_tol ? 1 : b_norm);

  const Real r_norm = sqrt<Real>(inner_prod(r, r, comm_));
  for (Long i = 0; i < N; i++) Q(0,i) = r[i] / r_norm;
  Vector<Real> beta(1); beta = r_norm;
  Vector<Real> sn, cs, h_k, q_k(N);

  Long k = 0;
  Real error = r_norm;
  for (; k < max_iter && error > abs_tol; k++) {
    if (verbose_ && !comm_.Rank()) printf("%3lld KSP Residual norm %.12e\n", (long long)k, (double)error);
    if (sn.Dim() <= k) ResizeVector(sn, (Long)((k+1)*ARRAY_RESIZE_FACTOR));
    if (cs.Dim() <= k) ResizeVector(cs, (Long)((k+1)*ARRAY_RESIZE_FACTOR));
    if (beta.Dim() <= k+1) ResizeVector(beta, (Long)((k+2)*ARRAY_RESIZE_FACTOR));
    if ( h_k.Dim() <= k+1) ResizeVector( h_k, (Long)((k+2)*ARRAY_RESIZE_FACTOR));

    arnoldi(h_k, q_k, A, k);
    apply_givens_rotation(h_k, cs[k], sn[k], cs, sn, k); // eliminate the last element in H ith row and update the rotation matrix
    for (Long i = 0; i < k+1; i++) H(k,i) = h_k[i];
    for (Long i = 0; i < N; i++) Q(k+1,i) = q_k[i];

    // update the residual vector
    beta[k+1] = -sn[k] * beta[k];
    beta[k]   = cs[k] * beta[k];
    error     = fabs(beta[k+1]);
  }
  if (verbose_ && !comm_.Rank()) printf("%3lld KSP Residual norm %.12e\n", (long long)k, (double)error);

  for (Long i = k-1; i >= 0; i--) { // beta <-- beta * inv(H); (through back substitution)
    beta[i] /= H(i,i);
    for (Long j = 0; j < i; j++) {
      beta[j] -= beta[i] * H(i,j);
    }
  }
  Vector<Real> x_(N); x_ = 0;
  for (Long i = 0; i < N; i++) { // x <-- beta * Q
    for (Long j = 0; j < k; j++) {
      x_[i] += beta[j] * Q(j,i);
    }
  }
  if (precond) precond->Apply(x_);
  (*x) += x_;

  if (solve_iter) (*solve_iter) = k;
  if (precond) std::cout<<"GMRES iter = "<<k<<", precond-size = "<<precond->Rank()<<'\n'; //////////////////////////////////////
  else std::cout<<"GMRES iter = "<<k<<'\n'; //////////////////////////////////////
  { ////////////////////////
    static Long fn_cnt = 0;
    fn_cnt += k;
    std::cout<<"GMRES fn evals = "<<fn_cnt<<'\n';
  }

  if (precond) {
    if (0) {
      Matrix<Real> Qt(N, k), U(k, N);
      for (Long i = 0; i < N; i++) { // apply givens rotations to Qt
        for (Long j = 0; j < k; j++) Qt[i][j] = Q_mat[j*N+i];
        for (Long j = 0; j < k-1; j++) { // apply givens rotations to Qt
          Real temp = cs[j] * Qt[i][j] + sn[j] * Qt[i][j+1];
          Qt[i][j+1] = -sn[j] * Qt[i][j] + cs[j] * Qt[i][j+1];
          Qt[i][j] = temp;
        }
      }
      for (Long i = 0; i < N; i++) {
        Qt[i][k-1] = cs[k-1] * Qt[i][k-1] + sn[k-1] * Q_mat[k*N+i];
      }

      Matrix<Real> Hinv(k,k);
      for (Long l = 0; l < k; l++) {
        for (Long i = 0; i < k; i++) Hinv[l][i] = 0;
        Hinv[l][l] = 1;
        for (Long i = l; i >= 0; i--) {
          Hinv[l][i] /= H(i,i);
          for (Long j = 0; j < i; j++) {
            Hinv[l][j] -= Hinv[l][i] * H(i,j);
          }
        }
      }
      Matrix<Real>::GEMM(U, Hinv, Matrix<Real>(k, N, Q_mat.begin(), false));
      U -= Qt.Transpose();

      precond->Append(Qt, U);
    } else {
      //if (precond->Size() == N) {
      //  const auto x = precond->Get_x();
      //  const auto Ax = precond->Get_Ax();
      //  for (Long k = 0; k < x.Dim()/N; k++) {
      //    precond_.RankOneUpdate(Vector<Real>(N,(Iterator<Real>)Ax.begin()+k*N,false), Vector<Real>(N,(Iterator<Real>)x.begin()+k*N,false));
      //  }
      //}
      //precond_.AugmentSubspace(*precond);
      //(*precond) = precond_;
      precond_.SetRHS(b);
      precond->AugmentSubspace(precond_);
    }
  }
}

template <class Real> inline void GMRES<Real>::operator()(Vector<Real>* x, const ParallelOp& A, const Vector<Real>& b, const Real tol, const Integer max_iter, const bool use_abs_tol, Long* solve_iter, KrylovPrecond<Real>* precond) {
  GenericGMRES(x, A, b, tol, max_iter, use_abs_tol, solve_iter, precond);
}

}  // end namespace

#ifdef SCTL_HAVE_PETSC

#include <petscksp.h>

namespace SCTL_NAMESPACE {

template <class Real> int GMRESMatVec(Mat M_, ::Vec x_, ::Vec Mx_) {
  PetscErrorCode ierr;

  PetscInt N, N_;
  VecGetLocalSize(x_, &N);
  VecGetLocalSize(Mx_, &N_);
  SCTL_ASSERT(N == N_);

  void* data = nullptr;
  MatShellGetContext(M_, &data);
  auto& M = dynamic_cast<const typename GMRES<Real>::ParallelOp&>(*(typename GMRES<Real>::ParallelOp*)data);

  const PetscScalar* x_ptr;
  ierr = VecGetArrayRead(x_, &x_ptr);
  CHKERRQ(ierr);

  Vector<Real> x(N);
  for (Long i = 0; i < N; i++) x[i] = (Real)x_ptr[i];
  Vector<Real> Mx(N);
  M(&Mx, x);

  PetscScalar* Mx_ptr;
  ierr = VecGetArray(Mx_, &Mx_ptr);
  CHKERRQ(ierr);

  for (long i = 0; i < N; i++) Mx_ptr[i] = Mx[i];
  ierr = VecRestoreArray(Mx_, &Mx_ptr);
  CHKERRQ(ierr);

  return 0;
}

PetscErrorCode MyKSPMonitor(KSP ksp, PetscInt n, PetscReal rnorm, void *dummy) {
  Comm* comm = (Comm*)dummy;
  if (!comm->Rank()) printf("%3lld KSP Residual norm %.12e\n", (long long)n, (double)rnorm);
  //PetscPrintf(PETSC_COMM_WORLD,"iteration %D KSP Residual norm %14.12e \n",n,rnorm);

  //PetscViewerAndFormat *vf;
  //PetscViewerAndFormatCreate(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_DEFAULT, &vf);
  //KSPMonitorResidual(ksp, n, rnorm, vf);
  //PetscViewerAndFormatDestroy(&vf);
  return 0;
}

template <class Real> inline void PETScGMRES(Vector<Real>* x, const typename GMRES<Real>::ParallelOp& A, const Vector<Real>& b, const Real tol, Integer max_iter, const bool use_abs_tol, const bool verbose_, const Comm& comm_, Long* solve_iter) {
  PetscInt N = b.Dim();
  if (max_iter < 0) { // set max_iter
    StaticArray<Long,2> NN{N,0};
    comm_.Allreduce(NN+0, NN+1, 1, Comm::CommOp::SUM);
    max_iter = NN[1];
  }
  const MPI_Comm comm = comm_.GetMPI_Comm();
  PetscErrorCode ierr;

  Mat PetscA;
  {  // Create Matrix. PetscA
    MatCreateShell(comm, N, N, PETSC_DETERMINE, PETSC_DETERMINE, (void*)&A, &PetscA);
    MatShellSetOperation(PetscA, MATOP_MULT, (void (*)(void))GMRESMatVec<Real>);
  }

  // Create linear solver context
  KSP ksp;
  ierr = KSPCreate(comm, &ksp);
  CHKERRABORT(comm, ierr);

  ::Vec Petsc_x, Petsc_b;
  {  // Create vectors
    VecCreateMPI(comm, N, PETSC_DETERMINE, &Petsc_b);
    VecCreateMPI(comm, N, PETSC_DETERMINE, &Petsc_x);

    PetscScalar* b_ptr;
    ierr = VecGetArray(Petsc_b, &b_ptr);
    CHKERRABORT(comm, ierr);
    for (long i = 0; i < N; i++) b_ptr[i] = b[i];
    ierr = VecRestoreArray(Petsc_b, &b_ptr);
    CHKERRABORT(comm, ierr);

    if (x->Dim() != N) {
      x->ReInit(N);
    } else {
      PetscScalar* x_ptr;
      ierr = VecGetArray(Petsc_x, &x_ptr);
      CHKERRABORT(comm, ierr);
      for (long i = 0; i < N; i++) x_ptr[i] = (*x)[i];
      ierr = VecRestoreArray(Petsc_x, &x_ptr);
      CHKERRABORT(comm, ierr);

      ierr = KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
    }
  }

  // Set operators. Here the matrix that defines the linear system
  // also serves as the preconditioning matrix.
  ierr = KSPSetOperators(ksp, PetscA, PetscA);
  CHKERRABORT(comm, ierr);

  // Set runtime options
  KSPSetType(ksp, KSPGMRES);
  KSPSetNormType(ksp, KSP_NORM_UNPRECONDITIONED);
  if (use_abs_tol) KSPSetTolerances(ksp, PETSC_DEFAULT, tol, PETSC_DEFAULT, max_iter);
  else KSPSetTolerances(ksp, tol, PETSC_DEFAULT, PETSC_DEFAULT, max_iter);
  KSPGMRESSetOrthogonalization(ksp, KSPGMRESModifiedGramSchmidtOrthogonalization);
  if (verbose_) KSPMonitorSet(ksp, MyKSPMonitor, comm, nullptr);
  KSPGMRESSetRestart(ksp, max_iter);
  ierr = KSPSetFromOptions(ksp);
  CHKERRABORT(comm, ierr);

  // -------------------------------------------------------------------
  // Solve the linear system: Ax=b
  // -------------------------------------------------------------------
  ierr = KSPSolve(ksp, Petsc_b, Petsc_x);
  CHKERRABORT(comm, ierr);

  // View info about the solver
  // KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD); CHKERRABORT(comm, ierr);

  // Iterations
  PetscInt its;
  ierr = KSPGetIterationNumber(ksp,&its); CHKERRABORT(comm, ierr);
  // ierr = PetscPrintf(PETSC_COMM_WORLD,"Iterations %D\n",its); CHKERRABORT(comm, ierr);

  {  // Set x
    const PetscScalar* x_ptr;
    ierr = VecGetArrayRead(Petsc_x, &x_ptr);
    CHKERRABORT(comm, ierr);

    for (long i = 0; i < N; i++) (*x)[i] = (Real)x_ptr[i];
  }

  ierr = KSPDestroy(&ksp);
  CHKERRABORT(comm, ierr);
  ierr = MatDestroy(&PetscA);
  CHKERRABORT(comm, ierr);
  ierr = VecDestroy(&Petsc_x);
  CHKERRABORT(comm, ierr);
  ierr = VecDestroy(&Petsc_b);
  CHKERRABORT(comm, ierr);

  if (solve_iter) (*solve_iter) = its;
}

template <> inline void GMRES<double>::operator()(Vector<double>* x, const ParallelOp& A, const Vector<double>& b, const double tol, const Integer max_iter, const bool use_abs_tol, Long* solve_iter) {
  PETScGMRES(x, A, b, tol, max_iter, use_abs_tol, verbose_, comm_, solve_iter);
}

template <> inline void GMRES<float>::operator()(Vector<float>* x, const ParallelOp& A, const Vector<float>& b, const float tol, const Integer max_iter, const bool use_abs_tol, Long* solve_iter) {
  PETScGMRES(x, A, b, tol, max_iter, use_abs_tol, verbose_, comm_, solve_iter);
}

}  // end namespace

#endif

#endif  //_SCTL_LIN_SOLVE_HPP_
