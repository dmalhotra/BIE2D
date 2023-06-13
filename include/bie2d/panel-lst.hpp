#ifndef BIE2D_PANEL_LST_HPP
#define BIE2D_PANEL_LST_HPP

#include <sctl.hpp>

namespace sctl {

  template <class Real, Integer Order, Integer digits, class Derived> class PanelLst {
    static constexpr Integer COORD_DIM = 2;

    public:

    static const Vector<Real>& Nds() {
      static const Vector<Real> nds = LegQuadRule<Real>::template nds<Order>();
      return nds;
    }
    static const Vector<Real>& Wts() {
      static const Vector<Real> wts = LegQuadRule<Real>::template wts<Order>();
      return wts;
    }

    PanelLst(const Vector<Real>& X = Vector<Real>()) {
      Init(X);
    }

    Long PanelCount() const {
      return Npanel;
    }

    void GetGeom(Vector<Real>* X = nullptr, Vector<Real>* Normal = nullptr, Vector<Real>* SurfWts = nullptr) const {
      if (X) *X = X_;
      if (Normal) *Normal = Normal_;
      if (SurfWts) {
        if (SurfWts->Dim() != Npanel*Order) SurfWts->ReInit(Npanel*Order);
        for (Long i = 0; i < Npanel; i++) {
          for (Long j = 0; j < Order; j++) {
            (*SurfWts)[i] = dS_[i*Order+j] * Wts()[j];
          }
        }
      }
    }

    void BoundaryIntegralDirect(Vector<Real>& I, const Vector<Real>& F) const {
      const Long dof = F.Dim() / (Npanel * Order);
      SCTL_ASSERT(F.Dim() == Npanel * Order * dof);

      if (I.Dim() != dof) I.ReInit(dof);
      for (Long k = 0; k < dof; k++) {
        Real I_ = 0;
        for (Long i = 0; i < Npanel; i++) {
          for (Long j = 0; j < Order; j++) {
            const Long idx = i*Order+j;
            I_ += F[idx*dof+k] * dS_[idx] * Wts()[j];
          }
        }
        I[k] = I_;
      }
    }

    template <class KerFn> void LayerPotential(Vector<Real>& U, const Vector<Real>& Xt, const Vector<Real>& F, const Real tol) const {
      static constexpr Integer DIM   = KerFn::CoordDim();
      static constexpr Integer KDIM0 = KerFn::SrcDim();
      static constexpr Integer KDIM1 = KerFn::TrgDim();
      SCTL_ASSERT(F.Dim() == Npanel * Order * KDIM0);
      static_assert(DIM == COORD_DIM, "Coordinate dimension mismatch.");

      const Long Nt = Xt.Dim() / COORD_DIM;
      if (U.Dim() != Nt * KDIM1) U.ReInit(Nt * KDIM1);

      #pragma omp parallel for schedule(static)
      for (Long t = 0; t < Nt; t++) {
        const Tensor<Real,false,COORD_DIM> Xt_((Iterator<Real>)Xt.begin()+t*COORD_DIM);
        Tensor<Real,false,KDIM1> UU(U.begin() + t*KDIM1);
        UU = (Real)0;

        for (Long i = 0; i < Npanel; i++) {
          Tensor<Real,true,Order,COORD_DIM> XX;
          for (Long j = 0; j < Order; j++) { // Set XX
            for (Integer k = 0; k < COORD_DIM; k++) {
              XX(j,k) = X_[(i*Order+j)*COORD_DIM+k] - Xt_(k);
            }
          }

          Tensor<Real,true,Order*KDIM0,KDIM1> MM;
          PanelKernelMatAdap<KerFn>(MM, XX, tol);

          const Tensor<Real,false,Order*KDIM0> FF((Iterator<Real>)F.begin() + i*Order*KDIM0);
          for (Integer k0 = 0; k0 < Order*KDIM0; k0++) {
            for (Integer k1 = 0; k1 < KDIM1; k1++) {
              UU(k1) += FF(k0) * MM(k0,k1);
            }
          }
        }
      }
    }

    template <class KerFn> void LayerPotentialMatrix(Matrix<Real>& M, const Vector<Real>& Xt, const Real tol) const {
      static constexpr Integer DIM   = KerFn::CoordDim();
      static constexpr Integer KDIM0 = KerFn::SrcDim();
      static constexpr Integer KDIM1 = KerFn::TrgDim();
      static_assert(DIM == COORD_DIM, "Coordinate dimension mismatch.");

      const Long Nt = Xt.Dim() / COORD_DIM;
      if (M.Dim(0) != Npanel*Order*KDIM0 || M.Dim(1) != Nt*KDIM1) {
        M.ReInit(Npanel*Order*KDIM0, Nt*KDIM1);
      }

      #pragma omp parallel for schedule(static)
      for (Long t = 0; t < Nt; t++) {
        const Tensor<Real,false,COORD_DIM> Xt_((Iterator<Real>)Xt.begin()+t*COORD_DIM);
        for (Long i = 0; i < Npanel; i++) {
          Tensor<Real,true,Order,COORD_DIM> XX;
          for (Long j = 0; j < Order; j++) { // Set XX
            for (Integer k = 0; k < COORD_DIM; k++) {
              XX(j,k) = X_[(i*Order+j)*COORD_DIM+k] - Xt_(k);
            }
          }

          Tensor<Real,true,Order*KDIM0,KDIM1> MM;
          PanelKernelMatAdap<KerFn>(MM, XX, tol);
          for (Long j = 0; j < Order*KDIM0; j++) {
            for (Long k = 0; k < KDIM1; k++) {
              M[i*Order*KDIM0+j][t*KDIM1+k] = MM(j, k);
            }
          }
        }
      }
    }

    void Init(const Vector<Real>& X) {
      Npanel = X.Dim() / (Order * COORD_DIM);
      SCTL_ASSERT(X.Dim() == Npanel * Order * COORD_DIM);
      if (dX_.Dim() != Npanel*Order*COORD_DIM) dX_.ReInit(Npanel*Order*COORD_DIM);
      if (Normal_.Dim() != Npanel*Order*COORD_DIM) Normal_.ReInit(Npanel*Order*COORD_DIM);
      if (dS_.Dim() != Npanel*Order) dS_.ReInit(Npanel*Order);

      X_ = X;
      for (Long i = 0; i < Npanel; i++) {
        const Vector<Real> X__(Order * COORD_DIM, X_.begin() + i*Order*COORD_DIM, false);
        Vector<Real> dX__(Order * COORD_DIM, dX_.begin() + i*Order*COORD_DIM, false);
        ComputeDerivative(dX__, X__);
      }
      for (Long i = 0; i < Npanel; i++) {
        for (Long j = 0; j < Order; j++) {
          const Real dx = dX_[(i*Order+j)*COORD_DIM+0];
          const Real dy = dX_[(i*Order+j)*COORD_DIM+1];
          const Real dr2 = dx*dx + dy*dy;
          SCTL_ASSERT(dr2 > 0);
          const Real dr = sqrt<Real>(dr2);
          const Real dr_inv = 1/dr;

          dS_[i*Order+j] = dr;
          Normal_[(i*Order+j)*COORD_DIM+0] = -dy * dr_inv;
          Normal_[(i*Order+j)*COORD_DIM+1] =  dx * dr_inv;
        }
      }
    }

    private:

    static Matrix<Real> InterpMat(const Vector<Real>& src_nds, const Vector<Real>& trg_nds) {
      Matrix<Real> M(src_nds.Dim(), trg_nds.Dim());
      Vector<Real> M_(M.Dim(0)*M.Dim(1), M.begin(), false);
      LagrangeInterp<Real>::Interpolate(M_, src_nds, trg_nds);
      return M;
    }

    static void ComputeDerivative(Vector<Real>& dY, const Vector<Real>& Y) {
      const auto DerivMat = []() {
        Vector<Real> V(Order); V = 0;
        Matrix<Real> M(Order, Order); M = 0;
        for (Long i = 0; i < Order; i++) {
          V[i] = 1;
          Vector<Real> dV(Order, M[i], false);
          LagrangeInterp<Real>::Derivative(dV, V, Nds());
          V[i] = 0;
        }
        return M;
      };
      static const Matrix<Real> Mt = DerivMat().Transpose();

      const Long dof = Y.Dim() / Order;
      SCTL_ASSERT(Y.Dim() == Order * dof);
      if (dY.Dim() != Order * dof) dY.ReInit(Order * dof);

      const Matrix<Real> Y_(Order, dof, (Iterator<Real>)Y.begin(), false);
      Matrix<Real> dY_(Order, dof, dY.begin(), false);
      Matrix<Real>::GEMM(dY_, Mt, Y_);
    }

    template <class KerFn, class Mat, class CoordVec> static void PanelKernelMat(Mat& M, const CoordVec& Xs) { // assume Xt = 0
      static constexpr Integer DIM   = KerFn::CoordDim();
      static constexpr Integer KDIM0 = KerFn::SrcDim();
      static constexpr Integer KDIM1 = KerFn::TrgDim();
      static_assert(Mat::template Dim<0>() == Order * KDIM0, "Output matrix dimension mismatch");
      static_assert(Mat::template Dim<1>() == KDIM1, "Output matrix dimension mismatch");
      static_assert(DIM == COORD_DIM, "Coordinate dimension mismatch.");

      static_assert(CoordVec::template Dim<0>() == Order, "Coordinate vector dimension mismatch");
      static_assert(CoordVec::template Dim<1>() == COORD_DIM, "Coordinate vector dimension mismatch");

      Tensor<Real,true,Order,2> dXs;
      { // Set dXs
        Vector<Real> dX(Order*COORD_DIM, dXs.begin(), false);
        const Vector<Real> X(Order*COORD_DIM, (Iterator<Real>)Xs.begin(), false);
        ComputeDerivative(dX, X);
      }

      for (Long i = 0; i < Order; i++) {
        const Real r[2] = {-Xs(i,0), -Xs(i,1)};
        const Real dx[2] = {dXs(i,0), dXs(i,1)};
        const Real da = sqrt<Real>(dx[0]*dx[0] + dx[1]*dx[1]);
        const Real da_inv = 1/da;
        const Real n[2] = {-dx[1]*da_inv, dx[0]*da_inv};

        Real M_[KDIM0][KDIM1];
        KerFn::template uKerMatrix<digits,Real>(M_, r, n, nullptr);
        for (Long k0 = 0; k0 < KDIM0; k0++) {
          for (Long k1 = 0; k1 < KDIM1; k1++) {
            M(i*KDIM0+k0, k1) = M_[k0][k1] * da * Wts()[i] * KerFn::template uKerScaleFactor<Real>();
          }
        }
      }
    }

    template <class KerFn, class Mat, class CoordVec> static void PanelKernelMatSing(Mat& M, const CoordVec& Xs, const Integer idx, const Real tol) { // assume Xt = 0
      static constexpr Integer DIM   = KerFn::CoordDim();
      static constexpr Integer KDIM0 = KerFn::SrcDim();
      static constexpr Integer KDIM1 = KerFn::TrgDim();
      static_assert(Mat::template Dim<0>() == Order * KDIM0, "Output matrix dimension mismatch");
      static_assert(Mat::template Dim<1>() == KDIM1, "Output matrix dimension mismatch");
      static_assert(DIM == COORD_DIM, "Coordinate dimension mismatch.");

      static_assert(CoordVec::template Dim<0>() == Order, "Coordinate vector dimension mismatch");
      static_assert(CoordVec::template Dim<1>() == COORD_DIM, "Coordinate vector dimension mismatch");

      Tensor<Real,true,Order,2> dXs;
      { // Set dXs
        Vector<Real> dX(Order*COORD_DIM, dXs.begin(), false);
        const Vector<Real> X(Order*COORD_DIM, (Iterator<Real>)Xs.begin(), false);
        ComputeDerivative(dX, X);
      }

      static const Vector<Vector<Real>> nds_wts = []() {
        Vector<Vector<Real>> nds_wts(2*Order);
        for (Long i = 0; i < Order; i++) {
          const Integer RefLevels = 2;
          const Integer LegQuadOrder = Order;

          auto DyadicQuad = [](Vector<Real>& nds, Vector<Real>& wts, const Integer LegQuadOrder, const Real s, const Integer levels, bool sort) {
            constexpr Integer LogQuadOrder = 16;
            //const auto& log_quad_nds = LogSingularityQuadRule<Real>(LogQuadOrder).first;
            //const auto& log_quad_wts = LogSingularityQuadRule<Real>(LogQuadOrder).second;
            const Real log_quad_nds[] = {
             (Real)2.0514967686824876939e-4L,
             (Real)1.0941116879700402194e-3L,
             (Real)4.3764467518801608777e-3L,
             (Real)1.3129579319567921241e-2L,
             (Real)3.1319422347023467220e-2L,
             (Real)6.4282044649094551893e-2L,
             (Real)1.1747685196991742595e-1L,
             (Real)1.9517904328816939977e-1L,
             (Real)2.9080420394081645443e-1L,
             (Real)4.0548336046512355400e-1L,
             (Real)5.4163934280979149964e-1L,
             (Real)6.8903327906975289198e-1L,
             (Real)8.1096672093024710801e-1L,
             (Real)8.9441573256030574193e-1L,
             (Real)9.7374799946956884221e-1L,
             (Real)9.9888924244762452447e-1L
            };
            const Real log_quad_wts[] = {
             (Real)6.4792916533086934623e-4L,
             (Real)1.2239247035761293717e-3L,
             (Real)5.9008349675174082927e-3L,
             (Real)1.2273728802285266544e-2L,
             (Real)2.5005188361723080279e-2L,
             (Real)4.1791379511932083369e-2L,
             (Real)6.5667543912496655405e-2L,
             (Real)8.8087941085870205392e-2L,
             (Real)1.0358776334798933346e-1L,
             (Real)1.2634697527479733119e-1L,
             (Real)1.4452810601255058206e-1L,
             (Real)1.4433961886609537948e-1L,
             (Real)9.2589501060143818582e-2L,
             (Real)8.7562563835839909214e-2L,
             (Real)5.7349021984199575159e-2L,
             (Real)3.0979791076523728382e-3L
            };
            const auto& leg_nds = LegendreQuadRule<Real>(LegQuadOrder).first;
            const auto& leg_wts = LegendreQuadRule<Real>(LegQuadOrder).second;

            Real len0 = std::min(pow<Real>(0.5,levels), std::min(s, (1-s)));
            Real len1 = std::min<Real>(s, 1-s);
            Real len2 = std::max<Real>(s, 1-s);

            for (Long i = 0; i < LogQuadOrder; i++) {
              nds.PushBack( len0*log_quad_nds[i]);
              nds.PushBack(-len0*log_quad_nds[i]);
              wts.PushBack(len0*log_quad_wts[i]);
              wts.PushBack(len0*log_quad_wts[i]);
            }

            for (Real start = len0; start < len1; start*=2) {
              Real step_ = std::min(start, len1-start);
              for (Long i = 0; i < leg_nds.Dim(); i++) {
                nds.PushBack( start + step_*leg_nds[i]);
                nds.PushBack(-start - step_*leg_nds[i]);
                wts.PushBack(step_*leg_wts[i]);
                wts.PushBack(step_*leg_wts[i]);
              }
            }

            for (Real start = len1; start < len2; start*=2) {
              Real step_ = std::min(start, len2-start);
              for (Long i = 0; i < leg_nds.Dim(); i++) {
                if (s + start + step_*leg_nds[i] <= 1.0) {
                  nds.PushBack( start + step_*leg_nds[i]);
                  wts.PushBack(step_*leg_wts[i]);
                }
                if (s - start - step_*leg_nds[i] >= 0.0) {
                  nds.PushBack(-start - step_*leg_nds[i]);
                  wts.PushBack(step_*leg_wts[i]);
                }
              }
            }

            if (!sort) return;
            Vector<Real> nds_(nds.Dim());
            Vector<Real> wts_(wts.Dim());
            Vector<std::pair<Real,Long>> sort_pair;
            for (Long i = 0; i < nds.Dim(); i++) {
              sort_pair.PushBack(std::pair<Real,Long>{nds[i], i});
            }
            std::sort(sort_pair.begin(), sort_pair.end());
            for (Long i = 0; i < nds.Dim(); i++) {
              const Long idx = sort_pair[i].second;
              nds_[i] = nds[idx];
              wts_[i] = wts[idx];
            }
            nds = nds_;
            wts = wts_;
          };

          auto& nds = nds_wts[2*i+0];
          auto& wts = nds_wts[2*i+1];
          DyadicQuad(nds, wts, LegQuadOrder, Nds()[i], RefLevels, false);
        }
        return nds_wts;
      }();
      static const Vector<Matrix<Real>> interp_mat = []() {
        Vector<Matrix<Real>> interp_mat(2*Order);
        for (Long i = 0; i < Order; i++) {
          interp_mat[2*i+0] = InterpMat(Nds()-Nds()[i], nds_wts[2*i+0]);
          interp_mat[2*i+1] = interp_mat[2*i+0].Transpose();
        }
        return interp_mat;
      }();

      //const Vector<Real>& nds = nds_wts[idx*2+0]; // in the interval [0,1]-Nds()[idx]
      const Vector<Real>& wts = nds_wts[idx*2+1];
      const Matrix<Real>& Minterp = interp_mat[idx*2+0];
      const Matrix<Real>& Minterp_t = interp_mat[idx*2+1];
      const Long Nnds = wts.Dim();

      Matrix<Real> Xs_(Nnds,COORD_DIM), dXs_(Nnds,COORD_DIM);
      Matrix<Real>::GEMM(Xs_, Minterp_t, Matrix<Real>(Order,COORD_DIM, (Iterator<Real>)Xs.begin(), false));
      Matrix<Real>::GEMM(dXs_, Minterp_t, Matrix<Real>(Order,COORD_DIM, (Iterator<Real>)dXs.begin(), false));

      Matrix<Real> MM(Nnds*KDIM0, KDIM1);
      for (Long i = 0; i < Nnds; i++) {
        const Real r[2] = {-Xs_[i][0], -Xs_[i][1]};
        const Real dx[2] = {dXs_[i][0], dXs_[i][1]};
        const Real da = sqrt<Real>(dx[0]*dx[0] + dx[1]*dx[1]);
        const Real da_inv = 1/da;
        const Real n[2] = {-dx[1]*da_inv, dx[0]*da_inv};

        Real M_[KDIM0][KDIM1];
        KerFn::template uKerMatrix<digits,Real>(M_, r, n, nullptr);
        for (Long k0 = 0; k0 < KDIM0; k0++) {
          for (Long k1 = 0; k1 < KDIM1; k1++) {
            MM(i*KDIM0+k0, k1) = M_[k0][k1] * da * wts[i] * KerFn::template uKerScaleFactor<Real>();
          }
        }
      }

      Matrix<Real> M_(Order, KDIM0*KDIM1, M.begin(), false);
      Matrix<Real>::GEMM(M_, Minterp, Matrix<Real>(Nnds, KDIM0*KDIM1, MM.begin(), false));
    }

    template <class KerFn, class Mat, class CoordVec> static void PanelKernelMatAdap(Mat& M, const CoordVec& Xs, const Real tol) { // assume Xt = 0
      static constexpr Integer DIM   = KerFn::CoordDim();
      static constexpr Integer KDIM0 = KerFn::SrcDim();
      static constexpr Integer KDIM1 = KerFn::TrgDim();
      static_assert(Mat::template Dim<0>() == Order * KDIM0, "Output matrix dimension mismatch");
      static_assert(Mat::template Dim<1>() == KDIM1, "Output matrix dimension mismatch");
      static_assert(DIM == COORD_DIM, "Coordinate dimension mismatch.");

      static_assert(CoordVec::template Dim<0>() == Order, "Coordinate vector dimension mismatch");
      static_assert(CoordVec::template Dim<1>() == COORD_DIM, "Coordinate vector dimension mismatch");

      static const Matrix<Real> Minterp1 = InterpMat(Nds(), Nds()*0.5);
      static const Matrix<Real> Minterp2 = InterpMat(Nds(), Nds()*0.5+0.5);
      static const Matrix<Real> Minterp1_t = Minterp1.Transpose();
      static const Matrix<Real> Minterp2_t = Minterp2.Transpose();

      Real min_dist;
      Integer min_idx;
      { // Set min_dist, min_idx
        Real min_dist2 = 0;
        for (Long i = 0; i < Order; i++) {
          const Real dist2 =  Xs(i,0)*Xs(i,0) + Xs(i,1)*Xs(i,1);
          if (i == 0 || dist2 < min_dist2) {
            min_dist2 = dist2;
            min_idx = i;
          }
        }
        min_dist = sqrt<Real>(min_dist2);
      }
      if (min_dist == 0) return PanelKernelMatSing<KerFn>(M, Xs, min_idx, tol);

      PanelKernelMat<KerFn>(M, Xs);
      if (tol <= 0) return;

      CoordVec Xs1, Xs2;
      Matrix<Real> Xs1_(Order,COORD_DIM, Xs1.begin(), false);
      Matrix<Real> Xs2_(Order,COORD_DIM, Xs2.begin(), false);
      const Matrix<Real> Xs_(Order,COORD_DIM, (Iterator<Real>)Xs.begin(), false);
      Matrix<Real>::GEMM(Xs1_, Minterp1_t, Xs_);
      Matrix<Real>::GEMM(Xs2_, Minterp2_t, Xs_);

      Mat MM;
      { // Set MM
        Mat Mchild1, Mchild2;
        PanelKernelMat<KerFn>(Mchild1, Xs1);
        PanelKernelMat<KerFn>(Mchild2, Xs2);
        const Matrix<Real> Mchild1_(Order, KDIM0*KDIM1, Mchild1.begin(), false);
        const Matrix<Real> Mchild2_(Order, KDIM0*KDIM1, Mchild2.begin(), false);
        Matrix<Real> MM_(Order, KDIM0*KDIM1, MM.begin(), false);
        MM_ = Minterp1 * Mchild1_ + Minterp2 * Mchild2_;
      };

      const Real rel_err = [&MM,&M]() {
        const Mat Merr = MM - M;
        const Matrix<Real> MM_(Order, KDIM0*KDIM1, (Iterator<Real>)MM.begin(), false);
        const Matrix<Real> Merr_(Order, KDIM0*KDIM1, (Iterator<Real>)Merr.begin(), false);
        Real max_err = 0, max_val = 0;
        for (Long k = 0; k < KDIM0*KDIM1; k++) {
          for (Long i = 0; i < Order; i++) {
            max_err = std::max<Real>(max_err, fabs(Merr_[i][k]));
            max_val = std::max<Real>(max_val, fabs(MM_[i][k]));
          }
        }
        return max_err/max_val;
      }();

      const Real panel_len = [&Xs]() {
        Tensor<Real,true,Order,2> dX;
        Vector<Real> dX_(Order*COORD_DIM, dX.begin(), false);
        ComputeDerivative(dX_, Vector<Real>(Order*COORD_DIM, (Iterator<Real>)Xs.begin(), false));

        Real sum = 0;
        for (Long i = 0; i < Order; i++) {
          const Real da = sqrt<Real>(dX(i,0)*dX(i,0) + dX(i,1)*dX(i,1));
          sum += da * Wts()[i];
        }
        return sum;
      }();

      if ((rel_err < tol || panel_len < min_dist/*TODO*/) && (min_dist > panel_len*(const_pi<Real>()/2/Order))) {
        M = MM;
        if (0) { // display error
          int depth = 0;
          Real t = panel_len/2/const_pi<Real>();
          while (t < sqrt<Real>(0.5)) {
            t *= 2;
            depth++;
            std::cout<<"  ";
          }
          std::cout<<depth<<' '<<rel_err<<'\n';
          if (depth > 150) {
            for (Long i = 0; i < Order; i++) {
              for (Long k = 0; k < COORD_DIM; k++) {
                std::cout<<Xs(i,k)<<' ';
              }
              std::cout<<'\n';
            }
            std::cout<<'\n';
            for (const auto& a : MM) std::cout<<a<<'\n';
            std::cout<<'\n';
            //for (const auto& a : Merr) std::cout<<a<<'\n';
            std::cout<<'\n';
            exit(0);
          }
        }
      } else {
        Mat M1, M2;
        PanelKernelMatAdap<KerFn>(M1, Xs1, tol);
        PanelKernelMatAdap<KerFn>(M2, Xs2, tol);
        const Matrix<Real> M1_(Order, KDIM0*KDIM1, M1.begin(), false);
        const Matrix<Real> M2_(Order, KDIM0*KDIM1, M2.begin(), false);
        Matrix<Real> M_(Order, KDIM0*KDIM1, M.begin(), false);
        M_ = Minterp1 * M1_ + Minterp2 * M2_;
      }
    }

    Long Npanel;
    Vector<Real> X_, dX_, Normal_; // Npanel * Order * COORD_DIM
    Vector<Real> dS_; // Npanel * Order

    friend Derived;
  };

}

#endif
