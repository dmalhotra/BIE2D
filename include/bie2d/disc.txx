namespace sctl {

  template <class Real, Integer Order, Integer digits> constexpr Integer Disc<Real,Order,digits>::CoordDim() {
    return COORD_DIM;
  }

  template <class Real, Integer Order, Integer digits> const Vector<Real>& Disc<Real,Order,digits>::PanelNds() {
    return PanelType::Nds();
  }

  template <class Real, Integer Order, Integer digits> const Vector<Real>& Disc<Real,Order,digits>::PanelWts() {
    return PanelType::Wts();
  }


  template <class Real, Integer Order, Integer digits> Disc<Real,Order,digits>::Disc() : radius_(0), coord_{0,0} {}

  template <class Real, Integer Order, Integer digits> Disc<Real,Order,digits>::Disc(const Real x, const Real y, const Real radius, const long Nunif) {
    Vector<Real> theta0, theta1;
    for (Long i = 0; i < Nunif; i++) {
      theta0.PushBack(2*const_pi<Real>()*(i+0)/Nunif);
      theta1.PushBack(2*const_pi<Real>()*(i+1)/Nunif);
    }
    Init(x, y, radius, theta0, theta1);
  }

  template <class Real, Integer Order, Integer digits> Disc<Real,Order,digits>::Disc(const Real x, const Real y, const Real radius, const Vector<Real>& theta0, const Vector<Real>& theta1) {
    Init(x, y, radius, theta0, theta1);
  }

  template <class Real, Integer Order, Integer digits> void Disc<Real,Order,digits>::Init(const Real x, const Real y, const Real radius, const Vector<Real>& theta0, const Vector<Real>& theta1) {
    const Long Npanel = theta0.Dim();
    SCTL_ASSERT((Long)theta1.Dim() == Npanel);

    coord_[0] = x;
    coord_[1] = y;
    radius_ = radius;
    theta0_ = theta0;
    theta1_ = theta1;

    const auto& nds = Disc::Nds();
    Vector<Real> X(Npanel * Order * COORD_DIM);
    for (Long i = 0; i < Npanel; i++) {
      for (Long j = 0; j < Order; j++) {
        const Real theta = theta0_[i] + (theta1_[i]-theta0_[i]) * nds[j];
        X[(i*Order+j)*2+0] = coord_[0] + radius * cos<Real>(theta);
        X[(i*Order+j)*2+1] = coord_[1] + radius * sin<Real>(theta);
      }
    }
    PanelType::Init(X);
  }



  template <class Real, Integer Order, Integer digits> Real Disc<Real,Order,digits>::Radius() const {
    return radius_;
  }

  template <class Real, Integer Order, Integer digits> Real Disc<Real,Order,digits>::Coord(int i) const {
    return coord_[i];
  }

  template <class Real, Integer Order, Integer digits> Long Disc<Real,Order,digits>::PanelCount() const {
    return theta0_.Dim();
  }

  template <class Real, Integer Order, Integer digits> Long Disc<Real,Order,digits>::NodeCount() const {
    return PanelCount() * Order;
  }

  template <class Real, Integer Order, Integer digits> std::pair<Real,Real> Disc<Real,Order,digits>::PanelRange(Long idx) const {
    return std::make_pair<Real,Real>(theta0_[idx], theta1_[idx]);
  }

  template <class Real, Integer Order, Integer digits> void Disc<Real,Order,digits>::GetGeom(Vector<Real>* X, Vector<Real>* Normal, Vector<Real>* SurfWts, Vector<Real>* theta) const {
    PanelType::GetGeom(X, Normal, SurfWts);
    if (theta) {
      const Long Npanel = PanelCount();
      if (theta->Dim() != Npanel * Order) theta->ReInit(Npanel * Order);

      for (Long i = 0; i < Npanel; i++) {
        for (Long j = 0; j < Order; j++) {
          (*theta)[i*Order+j] = theta0_[i] + (theta1_[i]-theta0_[i]) * PanelType::Nds()[j];
        }
      }
    }
  }



  template <class Real, Integer Order, Integer digits> void Disc<Real,Order,digits>::BoundaryIntegralDirect(Vector<Real>& I, const Vector<Real>& F) const {
    PanelType::BoundaryIntegralDirect(I, F);
  }

  template <class Real, Integer Order, Integer digits> template <class KerFn> void Disc<Real,Order,digits>::LayerPotential(Vector<Real>& U, const Vector<Real>& Xt, const Vector<Real>& F, const Real tol) const {
    PanelType::template LayerPotential<KerFn>(U, Xt, F, tol);
  }

  template <class Real, Integer Order, Integer digits> template <class KerFn> void Disc<Real,Order,digits>::LayerPotentialMatrix(Matrix<Real>& M, const Vector<Real>& Xt, const Real tol) const {
    PanelType::template LayerPotentialMatrix<KerFn>(M, Xt, tol);
  }



  template <class Disc> Long NodeCount(const Vector<Disc>& disc_lst) {
    Long N = 0;
    for (const auto& disc : disc_lst) N += disc.NodeCount();
    return N;
  }

  template <class Real, class Disc> void GetGeom(const Vector<Disc>& disc_lst, Vector<Real>* X, Vector<Real>* Normal, Vector<Real>* SurfWts) {
    const Long N = NodeCount(disc_lst);
    if (X && X->Dim() != N*Disc::CoordDim()) X->ReInit(N*Disc::CoordDim());
    if (Normal && Normal->Dim() != N*Disc::CoordDim()) Normal->ReInit(N*Disc::CoordDim());
    if (SurfWts && SurfWts->Dim() != N) SurfWts->ReInit(N);

    Long offset = 0;
    for (Long i = 0; i < disc_lst.Dim(); i++) {
      const Long N_ = disc_lst[i].NodeCount();
      Vector<Real> X_((X ? N_*Disc::CoordDim() : 0), (X ? X->begin() + offset*Disc::CoordDim() : NullIterator<Real>()), false);
      Vector<Real> Normal_((Normal ? N_*Disc::CoordDim() : 0), (Normal ? Normal->begin() + offset*Disc::CoordDim() : NullIterator<Real>()), false);
      Vector<Real> SurfWts_((SurfWts ? N_ : 0), (Normal ? Normal->begin() + offset : NullIterator<Real>()), false);
      disc_lst[i].GetGeom((X?&X_:nullptr), (Normal?&Normal_:nullptr), (SurfWts?&SurfWts_:nullptr));
      offset += N_;
    }
  }

  template <class Real, class Disc> void BoundaryIntegralDirect(Vector<Real>& I, const Vector<Disc>& disc_lst, const Vector<Real>& F) {
    const Long N = NodeCount(disc_lst);
    const Long dof = F.Dim() / N;
    SCTL_ASSERT(F.Dim() == N * dof);

    if (I.Dim() != dof) I.ReInit(dof);
    Vector<Real> I_(dof);
    I_ = 0;
    I = 0;
    Long offset = 0;
    for (const auto& disc : disc_lst) {
      const Long N_ = disc.NodeCount() * dof;
      const Vector<Real> F_(N_, (Iterator<Real>)F.begin() + offset, false);
      disc.BoundaryIntegralDirect(I_, F);
      I += I_;
      offset += N;
    }
  }

  template <class KerFn, class Real, class Disc> void LayerPotentialMatrix(Matrix<Real>& M, const Vector<Disc>& disc_lst, const Vector<Real>& Xt, const Real tol) {
    const Long Ns = NodeCount(disc_lst);
    const Long Nt = Xt.Dim() / KerFn::CoordDim();
    SCTL_ASSERT(Xt.Dim() == Nt * KerFn::CoordDim());

    if (M.Dim(0) != Ns * KerFn::SrcDim() || M.Dim(1) != Nt * KerFn::TrgDim()) {
      M.ReInit(Ns * KerFn::SrcDim(), Nt * KerFn::TrgDim());
    }

    Long src_offset = 0;
    for (const auto& disc : disc_lst) {
      const Long Ns_ = disc.NodeCount() * KerFn::SrcDim();
      Matrix<Real> M_(Ns_, Nt * KerFn::TrgDim(), M[src_offset], false);
      disc.template LayerPotentialMatrix<KerFn>(M_, Xt, tol);
      src_offset += Ns_;
    }
  }

  template <class KerFn, class Real, class Disc> void LayerPotential(Vector<Real>& U, const Vector<Disc>& disc_lst, const Vector<Real>& Xt, const Vector<Real>& F, const Real tol) {
    const Long Ns = NodeCount(disc_lst);
    const Long Nt = Xt.Dim() / KerFn::CoordDim();
    SCTL_ASSERT(Xt.Dim() == Nt * KerFn::CoordDim());
    SCTL_ASSERT(F.Dim() == Ns * KerFn::SrcDim());

    if (U.Dim() != Nt * KerFn::TrgDim()) U.ReInit(Nt * KerFn::TrgDim());
    Vector<Real> U_(Nt * KerFn::TrgDim());
    U_ = 0;
    U = 0;

    Long src_offset = 0;
    for (const auto& disc : disc_lst) {
      const Long Ns_ = disc.NodeCount() * KerFn::SrcDim();
      const Vector<Real> F_(Ns_, (Iterator<Real>)F.begin() + src_offset, false);
      disc.template LayerPotential<KerFn>(U_, Xt, F_, tol);
      U += U_;
      src_offset += Ns_;
    }

    //Matrix<Real> M;
    //LayerPotentialMatrix<Real,KerFn>(M, disc_lst, Xt, tol);
    //Matrix<Real> U_ = Matrix<Real>(1, F.Dim(), (Iterator<Real>)F.begin(), false) * M;
    //U = Vector<Real>(U_.Dim(1), U_.begin(), false);
  }

}
