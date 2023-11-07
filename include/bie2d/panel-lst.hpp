#ifndef BIE2D_PANEL_LST_HPP
#define BIE2D_PANEL_LST_HPP

#include <sctl.hpp>

namespace sctl {

  template <class Real, Integer Order> class PanelLst : public ElementListBase<Real> {
    static constexpr Integer COORD_DIM = 2;

    public:

    static constexpr Integer CoordDim() {
      return COORD_DIM;
    }
    static constexpr Integer ElemOrder() {
      return Order;
    }

    static const Vector<Real>& PanelNds() {
      static const Vector<Real> nds = LegQuadRule<Real>::template nds<Order>();
      return nds;
    }
    static const Vector<Real>& PanelWts() {
      static const Vector<Real> wts = LegQuadRule<Real>::template wts<Order>();
      return wts;
    }

    virtual ~PanelLst() {}

    PanelLst(const Vector<Real>& X = Vector<Real>()) {
      Init(X);
    }

    void Init(const Vector<Real>& X) {
      Npanel = X.Dim() / (Order * COORD_DIM);
      SCTL_ASSERT(X.Dim() == Npanel * Order * COORD_DIM);
      if (dX_.Dim() != Npanel*Order*COORD_DIM) dX_.ReInit(Npanel*Order*COORD_DIM);
      if (Normal_.Dim() != Npanel*Order*COORD_DIM) Normal_.ReInit(Npanel*Order*COORD_DIM);
      if (SurfWts_.Dim() != Npanel*Order) SurfWts_.ReInit(Npanel*Order);

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

          SurfWts_[i*Order+j] = dr * PanelWts()[j];
          Normal_[(i*Order+j)*COORD_DIM+0] =  dy * dr_inv;
          Normal_[(i*Order+j)*COORD_DIM+1] = -dx * dr_inv;
        }
      }
    }

    /**
     * Return the number of elements in the list.
     */
    virtual Long Size() const {
      return Npanel;
    }

    const Vector<Real>& SurfCoord() const {
      return X_;
    }
    const Vector<Real>& SurfNormal() const {
      return Normal_;
    }
    const Vector<Real>& SurfWts() const {
      return SurfWts_;
    }

    void BoundaryIntegralDirect(Vector<Real>& I, const Vector<Real>& F, const Long panel_start = 0, const Long panel_end = -1) const {
      const Long panel_end_ = (panel_end == -1 ? Npanel : panel_end);
      SCTL_ASSERT(0 <= panel_start && panel_start <= panel_end_ && panel_end_ <= Npanel);
      const Long Npanel_ = panel_end_ - panel_start;

      const Long dof = F.Dim() / (Npanel_ * Order);
      SCTL_ASSERT(F.Dim() == Npanel_ * Order * dof);

      if (I.Dim() != dof) I.ReInit(dof);
      for (Long k = 0; k < dof; k++) {
        Real I_ = 0;
        for (Long i = 0; i < Npanel; i++) {
          for (Long j = 0; j < Order; j++) {
            const Long idx = i*Order+j;
            I_ += F[idx*dof+k] * SurfWts_[idx];
          }
        }
        I[k] = I_;
      }
    }

    template <class KerFn, Integer digits> void LayerPotentialMatrix(Matrix<Real>& M, const Vector<Real>& Xt, const Real tol, const Long panel_start = 0, const Long panel_end = -1) const {
      static constexpr Integer DIM   = KerFn::CoordDim();
      static constexpr Integer KDIM0 = KerFn::SrcDim();
      static constexpr Integer KDIM1 = KerFn::TrgDim();
      static_assert(DIM == COORD_DIM, "Coordinate dimension mismatch.");

      const Long Nt = Xt.Dim() / COORD_DIM;
      const Long panel_end_ = (panel_end == -1 ? Npanel : panel_end);
      SCTL_ASSERT(0 <= panel_start && panel_start <= panel_end_ && panel_end_ <= Npanel);
      const Long Npanel_ = panel_end_ - panel_start;

      if (M.Dim(0) != Npanel_*Order*KDIM0 || M.Dim(1) != Nt*KDIM1) {
        M.ReInit(Npanel_*Order*KDIM0, Nt*KDIM1);
      }

      #pragma omp parallel for schedule(static) collapse(2)
      for (Long t = 0; t < Nt; t++) {
        for (Long i = panel_start; i < panel_end_; i++) {
          const Tensor<Real,false,COORD_DIM> Xt_((Iterator<Real>)Xt.begin()+t*COORD_DIM);
          Tensor<Real,true,Order,COORD_DIM> XX;
          for (Long j = 0; j < Order; j++) { // Set XX
            for (Integer k = 0; k < COORD_DIM; k++) {
              XX(j,k) = X_[(i*Order+j)*COORD_DIM+k] - Xt_(k);
            }
          }

          Tensor<Real,true,Order*KDIM0,KDIM1> MM;
          PanelKernelMatAdap<KerFn,digits>(MM, XX, tol);
          for (Long j = 0; j < Order*KDIM0; j++) {
            for (Long k = 0; k < KDIM1; k++) {
              M[(i-panel_start)*Order*KDIM0+j][t*KDIM1+k] = MM(j, k);
            }
          }
        }
      }
    }


    /**
     * Returns the position and normals of the surface nodal points for each
     * element.
     *
     * @param[out] X the position of the node points in array-of-struct
     * order: {x_1, y_1, z_1, x_2, ..., x_n, y_n, z_n}
     *
     * @param[out] Xn the normal vectors of the node points in
     * array-of-struct order: {nx_1, ny_1, nz_1, nx_2, ..., nx_n, ny_n, nz_n}
     *
     * @param[out] element_wise_node_cnt the number of node points
     * belonging to each element.
     */
    virtual void GetNodeCoord(Vector<Real>* X, Vector<Real>* Xn, Vector<Long>* element_wise_node_cnt) const {
      if (X) (*X) = X_;
      if (Xn) (*Xn) = Normal_;
      if (element_wise_node_cnt) {
        if (element_wise_node_cnt->Dim() != Size()) element_wise_node_cnt->ReInit(Size());
        for (Long i = 0; i < Size(); i++) element_wise_node_cnt[0][i] = Order;
      }
    }

    /**
     * Given an accuracy tolerance, returns the quadrature node positions,
     * the normals at the nodes, the weights and the cut-off distance from
     * the nodes for computing the far-field potential from the surface (at
     * target points beyond the cut-off distance).
     *
     * @param[out] X the position of the quadrature node points in
     * array-of-struct order: {x_1, y_1, z_1, x_2, ..., x_n, y_n, z_n}
     *
     * @param[out] Xn the normal vectors at the quadrature node points in
     * array-of-struct order: {nx_1, ny_1, nz_1, nx_2, ..., nx_n, ny_n, nz_n}
     *
     * @param[out] wts the weights corresponding to each quadrature node.
     *
     * @param[out] dist_far the cut-off distance from each quadrature node
     * such that quadrature rule will be accurate to the specified tolerance
     * for target points further away than this distance.
     *
     * @param[out] element_wise_node_cnt the number of quadrature node
     * points belonging to each element.
     *
     * @param[in] tol the accuracy tolerance.
     */
    virtual void GetFarFieldNodes(Vector<Real>& X, Vector<Real>& Xn, Vector<Real>& wts, Vector<Real>& dist_far, Vector<Long>& element_wise_node_cnt, const Real tol) const {
      X = X_;
      Xn = Normal_;
      wts = SurfWts_;
      if (element_wise_node_cnt.Dim() != Size()) element_wise_node_cnt.ReInit(Size());
      for (Long i = 0; i < Size(); i++) element_wise_node_cnt[i] = Order;

      Vector<Real> dist_far_gauss(Order);
      for (Long i = 0; i < Order; i++) { // Set dist_far_gauss
        const Real rho=pow<Real>((64/(15*tol)), (1/(Real)(2*Order)));
        const Real a = (rho-1/rho)/4;
        const Real b = (rho+1/rho)/4;
        const Real c = 0.5;

        dist_far_gauss[i] = b - fabs(PanelNds()[i]-0.5);
        const Real cos_t = b * (PanelNds()[i]-0.5) / (c*c);
        if (fabs(cos_t) <= 1) dist_far_gauss[i] = a * sqrt<Real>(1 + ((a*a)/(b*b)-1) * cos_t*cos_t);
      }

      if (dist_far.Dim() != Size()*Order) dist_far.ReInit(Size()*Order);
      for (Long i = 0 ; i < Size(); i++) {
        for (Long j = 0; j < Order; j++) {
          const Long node_idx = i * Order + j;
          Real dxds = sqrt<Real>(dX_[node_idx*COORD_DIM+0]*dX_[node_idx*COORD_DIM+0] + dX_[node_idx*COORD_DIM+1]*dX_[node_idx*COORD_DIM+1]);
          dist_far[node_idx] = dist_far_gauss[j] * dxds*2;
        }
      }
    }

    /**
     * Interpolates the density from surface node points to far-field
     * quadrature node points.
     *
     * @param[out] Fout the interpolated density at far-field quadrature
     * nodes in array-of-struct order.
     *
     * @param[in] Fin the input density at surface node points in
     * array-of-struct order.
     */
    virtual void GetFarFieldDensity(Vector<Real>& Fout, const Vector<Real>& Fin) const {
      Fout = Fin;
    }

    /**
     * Apply the transpose of the GetFarFieldDensity() operator for a single
     * element applied to the column-vectors of Min and the result is
     * returned in Mout.
     *
     * @param[out] Mout the output matrix where the column-vectors are the
     * result of the application of the transpose operator.
     *
     * @param[in] Min the input matrix whose column-vectors are
     * multiplied by the transpose operator.
     *
     * @param[in] elem_idx the index of the element.
     */
    virtual void FarFieldDensityOperatorTranspose(Matrix<Real>& Mout, const Matrix<Real>& Min, const Long elem_idx) const {
      Mout = Min;
    }

    /**
     * Compute self-interaction operator for each element.
     *
     * @param[out] M_lst the vector of all self-interaction matrices
     * (in row-major format).
     *
     * @param[in] ker the kernel object.
     *
     * @param[in] tol the accuracy tolerance.
     *
     * @param[in] trg_dot_prod whether to compute dot product of the potential with the target-normal vector.
     *
     * @param[in] self pointer to element-list object.
     */
    template <class Kernel> static void SelfInterac(Vector<Matrix<Real>>& M_lst, const Kernel& ker, Real tol, bool trg_dot_prod, const ElementListBase<Real>* self) {
      const auto& panel_lst = *dynamic_cast<const PanelLst*>(self);
      const Long Npanel = panel_lst.Size();
      SCTL_ASSERT(!trg_dot_prod);

      if (M_lst.Dim() != Npanel) M_lst.ReInit(Npanel);
      if (tol > 1e-3) {
        #pragma omp parallel for schedule(static)
        for (Long i = 0; i < Npanel; i++) {
          const Vector<Real> Xt(Order*COORD_DIM, (Iterator<Real>)panel_lst.X_.begin() + i*Order*COORD_DIM, false);
          panel_lst.template LayerPotentialMatrix<Kernel,3>(M_lst[i], Xt, tol, i, i+1);
        }
      } else if (tol > 1e-6) {
        #pragma omp parallel for schedule(static)
        for (Long i = 0; i < Npanel; i++) {
          const Vector<Real> Xt(Order*COORD_DIM, (Iterator<Real>)panel_lst.X_.begin() + i*Order*COORD_DIM, false);
          panel_lst.template LayerPotentialMatrix<Kernel,6>(M_lst[i], Xt, tol, i, i+1);
        }
      } else if (tol > 1e-9) {
        #pragma omp parallel for schedule(static)
        for (Long i = 0; i < Npanel; i++) {
          const Vector<Real> Xt(Order*COORD_DIM, (Iterator<Real>)panel_lst.X_.begin() + i*Order*COORD_DIM, false);
          panel_lst.template LayerPotentialMatrix<Kernel,9>(M_lst[i], Xt, tol, i, i+1);
        }
      } else if (tol > 1e-12) {
        #pragma omp parallel for schedule(static)
        for (Long i = 0; i < Npanel; i++) {
          const Vector<Real> Xt(Order*COORD_DIM, (Iterator<Real>)panel_lst.X_.begin() + i*Order*COORD_DIM, false);
          panel_lst.template LayerPotentialMatrix<Kernel,12>(M_lst[i], Xt, tol, i, i+1);
        }
      } else if (tol > 1e-15) {
        #pragma omp parallel for schedule(static)
        for (Long i = 0; i < Npanel; i++) {
          const Vector<Real> Xt(Order*COORD_DIM, (Iterator<Real>)panel_lst.X_.begin() + i*Order*COORD_DIM, false);
          panel_lst.template LayerPotentialMatrix<Kernel,15>(M_lst[i], Xt, tol, i, i+1);
        }
      } else {
        #pragma omp parallel for schedule(static)
        for (Long i = 0; i < Npanel; i++) {
          const Vector<Real> Xt(Order*COORD_DIM, (Iterator<Real>)panel_lst.X_.begin() + i*Order*COORD_DIM, false);
          panel_lst.template LayerPotentialMatrix<Kernel,-1>(M_lst[i], Xt, tol, i, i+1);
        }
      }
    }

    /**
     * Compute near-interaction operator for a given element-idx and each target.
     *
     * @param[out] M the near-interaction matrix (in row-major format).
     *
     * @param[in] Xt the position of the target points in array-of-structure
     * order: {x_1, y_1, z_1, x_2, ..., x_n, y_n, z_n}
     *
     * @param[in] normal_trg the normal at the target points in array-of-structure
     * order: {nx_1, ny_1, nz_1, nx_2, ..., nx_n, ny_n, nz_n}
     *
     * @param[in] ker the kernel object.
     *
     * @param[in] tol the accuracy tolerance.
     *
     * @param[in] elem_idx the index of the source element.
     *
     * @param[in] self pointer to element-list object.
     */
    template <class Kernel> static void NearInterac(Matrix<Real>& M, const Vector<Real>& Xt, const Vector<Real>& normal_trg, const Kernel& ker, Real tol, const Long elem_idx, const ElementListBase<Real>* self) {
      const auto& panel_lst = *dynamic_cast<const PanelLst*>(self);
      if      (tol > 1e-03) panel_lst.template LayerPotentialMatrix<Kernel, 3>(M, Xt, tol, elem_idx, elem_idx+1);
      else if (tol > 1e-06) panel_lst.template LayerPotentialMatrix<Kernel, 6>(M, Xt, tol, elem_idx, elem_idx+1);
      else if (tol > 1e-09) panel_lst.template LayerPotentialMatrix<Kernel, 9>(M, Xt, tol, elem_idx, elem_idx+1);
      else if (tol > 1e-12) panel_lst.template LayerPotentialMatrix<Kernel,12>(M, Xt, tol, elem_idx, elem_idx+1);
      else if (tol > 1e-15) panel_lst.template LayerPotentialMatrix<Kernel,15>(M, Xt, tol, elem_idx, elem_idx+1);
      else                  panel_lst.template LayerPotentialMatrix<Kernel,-1>(M, Xt, tol, elem_idx, elem_idx+1);
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
          LagrangeInterp<Real>::Derivative(dV, V, PanelNds());
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

    template <class KerFn, Integer digits, class Mat, class CoordVec> static void PanelKernelMat(Mat& M, const CoordVec& Xs) { // assume Xt = 0
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
        using VecType = Vec<Real,1>; // TODO: vectorize
        const VecType r[2] = {-Xs(i,0), -Xs(i,1)};
        const Real dx[2] = {dXs(i,0), dXs(i,1)};
        const Real da = sqrt<Real>(dx[0]*dx[0] + dx[1]*dx[1]);
        const Real da_inv = 1/da;
        const VecType n[2] = {dx[1]*da_inv, -dx[0]*da_inv};

        VecType M_[KDIM0][KDIM1];
        KerFn::template uKerMatrix<digits>(M_, r, n, nullptr);
        for (Long k0 = 0; k0 < KDIM0; k0++) {
          for (Long k1 = 0; k1 < KDIM1; k1++) {
            M(i*KDIM0+k0, k1) = M_[k0][k1][0] * da * PanelWts()[i] * KerFn::template uKerScaleFactor<Real>();
          }
        }
      }
    }

    template <class KerFn, Integer digits, class Mat, class CoordVec> static void PanelKernelMatSing(Mat& M, const CoordVec& Xs, const Integer idx, const Real tol) { // assume Xt = 0
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
          const Integer RefLevels = 4; // TODO: this should depend on tol or digits
          constexpr Integer LegQuadOrder = Order;

          auto DyadicQuad = [](Vector<Real>& nds, Vector<Real>& wts, const Integer LegQuadOrder, const Real s, const Integer levels, bool sort) {
            constexpr Integer LogQuadOrder = 34;
            //const auto& log_quad_nds = LogSingularityQuadRule<Real>(LogQuadOrder).first;
            //const auto& log_quad_wts = LogSingularityQuadRule<Real>(LogQuadOrder).second;
            const Real log_quad_nds34[] = {
              atoreal<Real>("0.99589374098049266697596309872721578E+00"),
              atoreal<Real>("0.97846047244658663853849870104083622E+00"),
              atoreal<Real>("0.94748824468677748770114451050673144E+00"),
              atoreal<Real>("0.90363906555115476011448272204195993E+00"),
              atoreal<Real>("0.84787483794349543231853699514258199E+00"),
              atoreal<Real>("0.78145018983210690100654660304172346E+00"),
              atoreal<Real>("0.70591315487419388110254799805491972E+00"),
              atoreal<Real>("0.62311128432926166862553990243506078E+00"),
              atoreal<Real>("0.53520397015870045727005445771481967E+00"),
              atoreal<Real>("0.44467879105384520085628317016161344E+00"),
              atoreal<Real>("0.35436312886956254571903857945617996E+00"),
              atoreal<Real>("0.26741263065759020885765178159167977E+00"),
              atoreal<Real>("0.18724892011755124045339493437106812E+00"),
              atoreal<Real>("0.11741781144977654722460870434161616E+00"),
              atoreal<Real>("0.61353718894321290901407286040266580E-01"),
              atoreal<Real>("0.22063169692420577829539354877107843E-01"),
              atoreal<Real>("0.17734943783838474259145370216303788E-02"),
              atoreal<Real>("0.10813203952530700264242441755928084E-01"),
              atoreal<Real>("0.40684269193361971594223779345243266E-01"),
              atoreal<Real>("0.88586268270837482102923525488591479E-01"),
              atoreal<Real>("0.15181767995905495788844366587415198E+00"),
              atoreal<Real>("0.22711966623286574096726799862750778E+00"),
              atoreal<Real>("0.31097640323428139741276424416602541E+00"),
              atoreal<Real>("0.39988896714993113174924233372761941E+00"),
              atoreal<Real>("0.49056239170215919019871577607868784E+00"),
              atoreal<Real>("0.58000439994711341496093553341378493E+00"),
              atoreal<Real>("0.66555898834864181650937948247918575E+00"),
              atoreal<Real>("0.74490417999439976027874952548594873E+00"),
              atoreal<Real>("0.81603709662961469516357114783969153E+00"),
              atoreal<Real>("0.87725951876587011658872227523155104E+00"),
              atoreal<Real>("0.92716888895641309301581849572948771E+00"),
              atoreal<Real>("0.96465486705047101931847993966027589E+00"),
              atoreal<Real>("0.98889899334050901265469149134829594E+00"),
              atoreal<Real>("0.99937397236539300005642832208544304E+00")
            };
            const Real log_quad_wts34[] = {
              atoreal<Real>("0.10527015989629867530381828777103869E-01")/2,
              atoreal<Real>("0.24287035808453484121080804847225055E-01")/2,
              atoreal<Real>("0.37546502188451407184665558948417978E-01")/2,
              atoreal<Real>("0.49991470446808202193490030283235742E-01")/2,
              atoreal<Real>("0.61327889301882068749377963656700249E-01")/2,
              atoreal<Real>("0.71263500034872754389323369886362529E-01")/2,
              atoreal<Real>("0.79502853810907992280367954629583594E-01")/2,
              atoreal<Real>("0.85741280412453331008573909706248727E-01")/2,
              atoreal<Real>("0.89658988802631161657424723727961569E-01")/2,
              atoreal<Real>("0.90920015591137141411591029340806841E-01")/2,
              atoreal<Real>("0.89184404723597223286776198841422607E-01")/2,
              atoreal<Real>("0.84143999118285015169742128316343244E-01")/2,
              atoreal<Real>("0.75588491730211133782319674920314119E-01")/2,
              atoreal<Real>("0.63495965946936024061306150863410387E-01")/2,
              atoreal<Real>("0.48124910777281996680445229495945482E-01")/2,
              atoreal<Real>("0.30088721774603493572321646506928115E-01")/2,
              atoreal<Real>("0.11117071714120635818675390141603161E-01")/2,
              atoreal<Real>("0.20043859834001163176226329635125052E-01")/2,
              atoreal<Real>("0.39278268957573892332450140419996636E-01")/2,
              atoreal<Real>("0.56072599051345180635203539912449369E-01")/2,
              atoreal<Real>("0.69841590807486564502560060870191172E-01")/2,
              atoreal<Real>("0.80171744627945581265012118302599718E-01")/2,
              atoreal<Real>("0.86955461948381084402388351201233661E-01")/2,
              atoreal<Real>("0.90318922436727576027577437184523573E-01")/2,
              atoreal<Real>("0.90528728962789270380776106619732628E-01")/2,
              atoreal<Real>("0.87912714053773823303840228108091299E-01")/2,
              atoreal<Real>("0.82809745729914426517278218512809409E-01")/2,
              atoreal<Real>("0.75547075477021100316989955292124314E-01")/2,
              atoreal<Real>("0.66435870360283354908047968627190078E-01")/2,
              atoreal<Real>("0.55775266165766421424636826203876866E-01")/2,
              atoreal<Real>("0.43858410981297792739645541453296435E-01")/2,
              atoreal<Real>("0.30977161597657852070245070674234555E-01")/2,
              atoreal<Real>("0.17422913722490842876658619894040446E-01")/2,
              atoreal<Real>("0.35395471132811402225998935613870556E-02")/2
            };
            const Real log_quad_nds16[] = {
             atoreal<Real>("2.051496768682487693930555255204858e-4"),
             atoreal<Real>("1.094111687970040219444686100974433e-3"),
             atoreal<Real>("4.376446751880160877778744403897733e-3"),
             atoreal<Real>("1.312957931956792124115555363331109e-2"),
             atoreal<Real>("3.131942234702346722017523520478278e-2"),
             atoreal<Real>("6.414075003315194736170246864401867e-2"),
             atoreal<Real>("1.167656676788392962581590955470186e-1"),
             atoreal<Real>("1.907815000663038947234049372880373e-1"),
             atoreal<Real>("2.888708401162808885014121988374203e-1"),
             atoreal<Real>("4.062500000000000000000000000000000e-1"),
             atoreal<Real>("5.425366345565433699292444290664887e-1"),
             atoreal<Real>("6.772989490147958863901600197889881e-1"),
             atoreal<Real>("8.022989490147958863901600197889882e-1"),
             atoreal<Real>("9.024104783559153001128950056622488e-1"),
             atoreal<Real>("9.687500000000000000000000000000000e-1"),
             atoreal<Real>("9.998611553059530655596495295904344e-1")
            };
            const Real log_quad_wts16[] = {
             atoreal<Real>("6.257250862796849297416842333651581e-4"),
             atoreal<Real>("1.331513954762220548078572444722675e-3"),
             atoreal<Real>("5.630698641374881532106925242735762e-3"),
             atoreal<Real>("1.278437652722473355188049324544743e-2"),
             atoreal<Real>("2.414877914532726019371510435272053e-2"),
             atoreal<Real>("4.282309343925836415830433385515808e-2"),
             atoreal<Real>("6.209738255512643155810737848260408e-2"),
             atoreal<Real>("8.716870408981970882387058469424124e-2"),
             atoreal<Real>("1.071938042328450559331758423068564e-1"),
             atoreal<Real>("1.288829470706290078509662292589008e-1"),
             atoreal<Real>("1.387941899047568600794727699221984e-1"),
             atoreal<Real>("1.304683663335963765297430962804360e-1"),
             atoreal<Real>("1.163052762532103408290320175564605e-1"),
             atoreal<Real>("8.250707132002965824282082732808748e-2"),
             atoreal<Real>("5.032617692455981053749882887388336e-2"),
             atoreal<Real>("8.911894521199604701485311922178854e-3")
            };
            const auto& log_quad_nds = (LogQuadOrder==16 ? log_quad_nds16 : log_quad_nds34);
            const auto& log_quad_wts = (LogQuadOrder==16 ? log_quad_wts16 : log_quad_wts34);
            static const auto& leg_nds = LegQuadRule<Real>::ComputeNds(Order);
            static const auto& leg_wts = LegQuadRule<Real>::ComputeWts(leg_nds);

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
          DyadicQuad(nds, wts, LegQuadOrder, PanelNds()[i], RefLevels, false);
        }
        return nds_wts;
      }();
      static const Vector<Matrix<Real>> interp_mat = []() {
        Vector<Matrix<Real>> interp_mat(2*Order);
        for (Long i = 0; i < Order; i++) {
          interp_mat[2*i+0] = InterpMat(PanelNds()-PanelNds()[i], nds_wts[2*i+0]);
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
        using VecType = Vec<Real,1>; // TODO: vectorize
        const VecType r[2] = {-Xs_[i][0], -Xs_[i][1]};
        const Real dx[2] = {dXs_[i][0], dXs_[i][1]};
        const Real da = sqrt<Real>(dx[0]*dx[0] + dx[1]*dx[1]);
        const Real da_inv = 1/da;
        const VecType n[2] = {dx[1]*da_inv, -dx[0]*da_inv};

        VecType M_[KDIM0][KDIM1];
        KerFn::template uKerMatrix<digits>(M_, r, n, nullptr);
        for (Long k0 = 0; k0 < KDIM0; k0++) {
          for (Long k1 = 0; k1 < KDIM1; k1++) {
            MM(i*KDIM0+k0, k1) = M_[k0][k1][0] * da * wts[i] * KerFn::template uKerScaleFactor<Real>();
          }
        }
      }

      Matrix<Real> M_(Order, KDIM0*KDIM1, M.begin(), false);
      Matrix<Real>::GEMM(M_, Minterp, Matrix<Real>(Nnds, KDIM0*KDIM1, MM.begin(), false));
    }

    template <class KerFn, Integer digits, class Mat, class CoordVec> static void PanelKernelMatAdap(Mat& M, const CoordVec& Xs, const Real tol) { // assume Xt = 0
      static constexpr Integer DIM   = KerFn::CoordDim();
      static constexpr Integer KDIM0 = KerFn::SrcDim();
      static constexpr Integer KDIM1 = KerFn::TrgDim();
      static_assert(Mat::template Dim<0>() == Order * KDIM0, "Output matrix dimension mismatch");
      static_assert(Mat::template Dim<1>() == KDIM1, "Output matrix dimension mismatch");
      static_assert(DIM == COORD_DIM, "Coordinate dimension mismatch.");

      static_assert(CoordVec::template Dim<0>() == Order, "Coordinate vector dimension mismatch");
      static_assert(CoordVec::template Dim<1>() == COORD_DIM, "Coordinate vector dimension mismatch");

      static const Matrix<Real> Minterp1 = InterpMat(PanelNds(), PanelNds()*0.5);
      static const Matrix<Real> Minterp2 = InterpMat(PanelNds(), PanelNds()*0.5+0.5);
      static const Matrix<Real> Minterp1_t = Minterp1.Transpose();
      static const Matrix<Real> Minterp2_t = Minterp2.Transpose();

      Real min_dist;
      Integer min_idx = 0;
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
      if (min_dist == 0) return PanelKernelMatSing<KerFn,digits>(M, Xs, min_idx, tol);

      PanelKernelMat<KerFn,digits>(M, Xs);
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
        PanelKernelMat<KerFn,digits>(Mchild1, Xs1);
        PanelKernelMat<KerFn,digits>(Mchild2, Xs2);
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
          sum += da * PanelWts()[i];
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
        PanelKernelMatAdap<KerFn,digits>(M1, Xs1, tol);
        PanelKernelMatAdap<KerFn,digits>(M2, Xs2, tol);
        const Matrix<Real> M1_(Order, KDIM0*KDIM1, M1.begin(), false);
        const Matrix<Real> M2_(Order, KDIM0*KDIM1, M2.begin(), false);
        Matrix<Real> M_(Order, KDIM0*KDIM1, M.begin(), false);
        M_ = Minterp1 * M1_ + Minterp2 * M2_;
      }
    }

    Long Npanel;
    Vector<Real> X_, dX_, Normal_; // Npanel * Order * COORD_DIM
    Vector<Real> SurfWts_; // Npanel * Order
  };

}

#endif
