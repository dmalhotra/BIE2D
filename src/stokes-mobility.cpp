#include <bie2d.hpp>
using namespace sctl;

constexpr Integer ElemOrder = 24;
constexpr Integer COORD_DIM = 2;
using RefReal = QuadReal; // reference solution precision
using Real = double;

template <Integer ElemOrder, class Real> Vector<Real> DiscMobilitySolve(const Vector<Real>& F, const Vector<Real>& X, const Real R, const ICIPType icip_type, const Real tol) {
  const Comm comm = Comm::Self();
  DiscMobility<Real,ElemOrder> disc_mobility(comm);

  Vector<Real> V;
  disc_mobility.Init(X, R, tol, icip_type);
  disc_mobility.Solve(V, F, Vector<Real>(), 5000);
  return V;
}

template <class Real> void TestCapacElast(const ICIPType icip_type, const Real eps, const Real gmres_tol, const Long gmres_iter) {
  const Comm comm = Comm::Self();
  const Real R = 1; //0.75;
  const Real tol = machine_eps<Real>();
  Profile::Enable(false);

  const Long Ndisc = 2;
  Vector<Real> X(Ndisc*COORD_DIM), V0(Ndisc);
  { // Set V0, X
    //X[0] = -(R+eps/2)*sqrt<Real>(0.5); X[1] = -(R+eps/2)*sqrt<Real>(0.5);
    //X[2] =  (R+eps/2)*sqrt<Real>(0.5); X[3] =  (R+eps/2)*sqrt<Real>(0.5);
    X[0] = -(R+eps/2); X[1] = 0;
    X[2] =  (R+eps/2); X[3] = 0;
    //X[4] =  0; X[5] = sqrt<Real>((Real)3)*(R+eps/2);

    V0[0] = -0.5;
    V0[1] =  0.5;
    //V0[2] =  0.5;
  }

  DiscCapacitance<Real,ElemOrder> capacitance(comm);
  DiscElastance<Real,ElemOrder> elastance(comm);
  capacitance.Init(X, R, tol, icip_type);
  elastance.Init(X, R, tol, icip_type);

  Vector<Real> Q, V;
  capacitance.Solve(Q, V0, gmres_tol, gmres_iter);
  elastance.Solve(V, Q, gmres_tol, gmres_iter);
  std::cout<<V0<<Q<<V<<'\n';

  std::cout<<std::setprecision(18)<<Q[0]<<' '<<Q[1]<<'\n';

  Real max_err = 0;
  for (const auto& x : V0-V) max_err = std::max<Real>(max_err, fabs(x));
  std::cout<<"Error = "<<max_err<<'\n';

  //if (eps == 1e-6) std::cout<<"Capacitance-Error = "<<std::max(fabs(atoreal<Real>("-3.14159278448947793e+3")-Q[0]), fabs(atoreal<Real>("3.14159278448947793e+3")-Q[1]))/std::max(fabs(Q[0]), fabs(Q[1]))<<'\n';
  //if (eps == 1e-8) std::cout<<"Capacitance-Error = "<<std::max(fabs(atoreal<Real>("-3.14159265489879153e+4")-Q[0]), fabs(atoreal<Real>("3.14159265489879153e+4")-Q[1]))/std::max(fabs(Q[0]), fabs(Q[1]))<<'\n';

  Profile::print(&comm);
}

template <class Real> Vector<Real> test_mobility(const Real R, const Real eps, const ICIPType icip_type, const Real tol, const Real gmres_tol, const Long gmres_iter) {
  const Long Ndisc = 2;
  Vector<Real> X(Ndisc*COORD_DIM), F(Ndisc*3);
  { // Set F, X
    X = 0;
    X[0] = -(R+eps/2); X[1] = 0;
    X[2] =  (R+eps/2); X[3] = 0;
    if (Ndisc > 2) {
      X[4] =  0;
      X[5] = sqrt<Real>((Real)3)*(R+eps/2);
    }
    //X[0] = -(R+eps/2)*sqrt<Real>(0.5); X[1] = -(R+eps/2)*sqrt<Real>(0.5);
    //X[2] =  (R+eps/2)*sqrt<Real>(0.5); X[3] =  (R+eps/2)*sqrt<Real>(0.5);

    F = 0;
    //for (auto& f : F) f = drand48();
    F[0] = 0.5; F[1] = 1; F[2] = 1.3;
    F[3] = 0; F[4] = 0; F[5] = 0;
  }

  const Comm comm = Comm::Self();
  static DiscMobility<Real,ElemOrder> disc_mobility(comm);
  disc_mobility.Init(X, R, tol, icip_type);

  Vector<Real> V;
  disc_mobility.Solve(V, F, Vector<Real>(), gmres_tol, gmres_iter, nullptr);
  return V;
}

int main(int argc, char** argv) {
  Comm::MPI_Init(&argc, &argv);

  if (0) {
    const Real eps = 1e-8;

    TestCapacElast<QuadReal>(ICIPType::Adaptive, eps, 1e-22, 1000);

    //TestCapacElast<RefReal>(ICIPType::Compress, eps, 1e-22, 100);
    TestCapacElast<Real>(ICIPType::Compress, eps, 1e-22, 100);

    //TestCapacElast<RefReal>(ICIPType::Precond, eps, 1e-22, 50);
    TestCapacElast<Real>(ICIPType::Precond, eps, 1e-22, 50);

    return 0;
  }

  if (1) {
    const Real R = 0.75;
    const Real eps = (argc > 1 ? atoreal<Real>(argv[1]) : 1e-10);
    const Real tol = machine_eps<Real>()*64;

    std::cout<<"R = "<<R<<'\n';
    std::cout<<"eps = "<<eps<<'\n';

    const auto real2quad = [](const Vector<Real>& v) {
      Vector<RefReal> w(v.Dim());
      for (Long i = 0; i < v.Dim(); i++) w[i] = (RefReal)v[i];
      return w;
    };

    Vector<Vector<RefReal>> V(6); V = 0;
    const RefReal eps_ = pow<RefReal>((RefReal)10, round<RefReal>(log<RefReal>((RefReal)eps)/log<RefReal>(10)));
    SCTL_ASSERT(fabs(eps-eps_) < eps*1e-2);

    if (0) { ////////////////////////////////////////////
      V[0] = test_mobility<RefReal>(R, eps_, ICIPType::Adaptive, tol*1e-14, -1, 0);
      std::cout<<"Ref solution: eps="<<std::scientific<<std::setprecision(30)<<eps_<<";  V=";
      for (const auto& x : V[0]) std::cout<<std::setprecision(30)<<x<<' ';
      std::cout<<'\n';
      return 0;
    } else {
      if (fabs(eps-1e-01) < 1e-14) V[0] = Vector<RefReal>{1.92788066614632465125918674710e-2,  4.16196198920043460125827532648e-2,  6.65495736471997266941115534710e-2,  1.71798873865923935706407216610e-2,  2.30903030912668191519172611802e-2, -2.08291139309286707701873040494e-2};
      if (fabs(eps-1e-02) < 1e-14) V[0] = Vector<RefReal>{1.88462329076690489163492509559e-2,  1.88008975926765081963129146025e-2,  2.99407877186829776817227517333e-2,  1.87663913705439837757306129844e-2,  2.52896727880993981379450718936e-2, -2.48597500284358668417946226466e-3};
      if (fabs(eps-1e-03) < 1e-14) V[0] = Vector<RefReal>{1.88659916463826227080134474713e-2,  9.81363828081942325063300386670e-3,  1.55574002657018265519147543497e-2,  1.88634146771913954981858489708e-2,  2.38188470464957562628921194485e-2,  5.11773755250713813251504707769e-3};
      if (fabs(eps-1e-04) < 1e-14) V[0] = Vector<RefReal>{1.88705901749754363694457929278e-2,  6.89557410621338912536327777439e-3,  1.08706369338561432145964605774e-2,  1.88705085154023081327821022511e-2,  2.32173288657722692769754060820e-2,  7.56333113106109907438099402962e-3};
      if (fabs(eps-1e-05) < 1e-14) V[0] = Vector<RefReal>{1.88711353258985043667073642691e-2,  5.96957708032943417167055557957e-3,  9.38090544223907833389021253642e-3,  1.88711327430615401538401750202e-2,  2.30171292511774916553574237539e-2,  8.33485352730129599633643778668e-3};
      if (fabs(eps-1e-06) < 1e-14) V[0] = Vector<RefReal>{1.88711925449578704039704386569e-2,  5.67657211302266347879909813976e-3,  8.90924899523724756775061213119e-3,  1.88711924632797033864568771740e-2,  2.29529423343814540411860673965e-2,  8.57845232570321691887415119056e-3};
      if (fabs(eps-1e-07) < 1e-14) V[0] = Vector<RefReal>{1.88711983523812078827890158901e-2,  5.58390246110060206837798267963e-3,  8.76004836962730420952627080560e-3,  1.88711983497983121073754511385e-2,  2.29325606355103788002196017535e-2,  8.65544108777823492459351606779e-3};
      if (fabs(eps-1e-08) < 1e-14) V[0] = Vector<RefReal>{1.88711989358277674480722239714e-2,  5.55459655362608874252962763063e-3,  8.71286221670679268158017592263e-3,  1.88711989357460890950735243178e-2,  2.29261070955936016140367724090e-2,  8.67978248364906175167312828745e-3};
      if (fabs(eps-1e-09) < 1e-14) V[0] = Vector<RefReal>{1.88711989942579372570688838565e-2,  5.54532909748470724378051682262e-3,  8.69794017418499597836043612128e-3,  1.88711989942553543607053510812e-2,  2.29240654828132799147912114029e-2,  8.68747944390978939040427574196e-3};
      if (fabs(eps-1e-10) < 1e-14) V[0] = Vector<RefReal>{1.88711990001036584088926238595e-2,  5.54239845922764475389764944147e-3,  8.69322136313803770515531582436e-3,  1.88711990001035767305377647000e-2,  2.29234197858583061568077350956e-2,  8.68991338976621305466029742905e-3};
      if (fabs(eps-1e-11) < 1e-14) V[0] = Vector<RefReal>{1.88711990006883160373171493308e-2,  5.54147170891268796666748105149e-3,  8.69172913938126028332107227316e-3,  1.88711990006883134544207785183e-2,  2.29232155903258862652755173466e-2,  8.69068306635165816471724222223e-3};
      if (fabs(eps-1e-12) < 1e-14) V[0] = Vector<RefReal>{1.88711990007467845043242528941e-2,  5.54117864461848595707447833165e-3,  8.69125725632834742583242901372e-3,  1.88711990007467844226459933015e-2,  2.29231510172062768807399638748e-2,  8.69092645899110041045027476552e-3};

      // Interpolation interval: 1 - 1e-14
      //
      // double:
      //   eps   Adaptive        32         64         96        128        160        192        224        256        288        320
      // 1e-01  7.558e-15  7.256e-5   2.790e-8  6.274e-11  9.801e-15  9.383e-15  8.758e-15  7.298e-15  7.715e-15  8.966e-15  8.132e-15
      // 1e-02  1.800e-13  1.845e-6   3.667e-7  6.391e-10  1.493e-12  1.575e-14  1.239e-14  1.008e-14  1.263e-14  1.529e-14  1.019e-14
      // 1e-03  4.430e-13  3.172e-4   4.114e-6   2.581e-9  2.651e-11  1.913e-13  8.309e-14  6.256e-14  7.814e-14  7.989e-14  7.304e-14
      // 1e-04  3.800e-11  2.207e-2   1.234e-7   2.992e-8  2.215e-10  5.684e-13  6.515e-14  1.068e-13  1.405e-13  1.329e-13  9.922e-14
      // 1e-05   8.950e-9  9.033e-2   1.246e-4   2.301e-8  5.899e-10  7.164e-12  1.741e-12  8.614e-13  1.266e-12  1.593e-13  1.603e-12
      // 1e-06   1.952e-7  3.503e-1   1.106e-3   2.227e-6   4.102e-9  5.542e-11  1.642e-11  1.291e-11  1.404e-11  1.431e-11  8.302e-12
      // 1e-07   4.301e-7  9.804e-1   6.764e-3   2.828e-6   7.126e-8  6.380e-10  2.987e-10  1.101e-10  3.618e-10  1.584e-10  7.350e-11
      // 1e-08   5.274e-8  1.000e+0   3.845e-2   1.445e-4   8.533e-7   7.374e-9   1.944e-9   1.080e-8   1.110e-8   2.996e-9   5.307e-9
      // 1e-09   2.347e-7  1.000e+0   1.313e-1   7.355e-4   7.806e-6   4.682e-8   3.016e-8   8.668e-8   1.164e-7   1.591e-7   3.303e-8
      // 1e-10   9.592e-8  1.000e+0   5.393e-1   1.334e-2   4.112e-5   1.592e-6   8.762e-7   6.882e-7   4.735e-7   4.092e-7   8.456e-7
      // 1e-11   4.418e-9  9.999e-1   1.019e+0   1.044e-1   7.043e-4   4.338e-6   4.546e-6   3.393e-6   1.498e-5   1.034e-5   4.081e-6
      // 1e-12  4.495e-10  9.997e-1   1.005e+0   4.582e-1   2.600e-4   1.060e-4   6.369e-5   3.279e-5   1.333e-4   7.902e-5   5.898e-5
      //
      // QuadReal:
      //   eps                   32         64         96        128        160        192        224        256        288        320
      // 1e-01             7.256e-5   2.790e-8  6.274e-11  7.530e-15  2.652e-16  4.639e-17  4.358e-17  4.362e-17  4.362e-17  4.362e-17
      // 1e-02             1.845e-6   3.667e-7  6.391e-10  1.478e-12  1.090e-14  9.594e-17  3.109e-17  3.109e-17  3.109e-17  3.109e-17
      // 1e-03             3.172e-4   4.114e-6   2.581e-9  2.646e-11  1.317e-13  3.777e-16  2.825e-17  2.874e-17  2.874e-17  2.874e-17
      // 1e-04             2.207e-2   1.234e-7   2.992e-8  2.215e-10  5.245e-13  2.541e-15  2.673e-17  2.668e-17  2.668e-17  2.668e-17
      // 1e-05             9.033e-2   1.246e-4   2.301e-8  5.899e-10  6.829e-12  3.660e-14  9.895e-17  6.490e-17  6.490e-17  6.490e-17
      // 1e-06             3.503e-1   1.106e-3   2.227e-6   4.104e-9  5.458e-11  2.707e-13  3.039e-16  2.184e-16  2.257e-16  2.244e-16
      // 1e-07             9.804e-1   6.764e-3   2.828e-6   7.128e-8  6.156e-10  2.198e-12  4.667e-15  5.288e-16  5.700e-16  5.587e-16
      // 1e-08             1.000e+0   3.845e-2   1.445e-4   8.535e-7   1.659e-9  1.029e-11  6.818e-14  8.883e-15  8.427e-15  8.436e-15
      // 1e-09             1.000e+0   1.313e-1   7.355e-4   7.807e-6   4.444e-8  5.989e-11  2.480e-13  4.462e-14  4.002e-14  3.978e-14
      // 1e-10             1.000e+0   5.393e-1   1.334e-2   4.113e-5   1.716e-7   1.474e-9  4.056e-12  1.511e-13  1.462e-13  1.439e-13
      // 1e-11             9.999e-1   1.019e+0   1.044e-1   7.044e-4   4.278e-6   1.381e-8  5.144e-11  9.219e-13  4.478e-13  4.648e-13
      // 1e-12             1.000e+0   1.005e+0   4.597e-1   2.749e-4   1.794e-5   1.079e-7  5.166e-10  5.876e-12  1.615e-12  1.522e-12


      // Interpolation interval: 1 - 1e-8
      //
      // double:
      //   eps   Adaptive        32         64         96        128        160
      // 1e-01  7.558e-15  4.474e-6  1.420e-10  9.592e-15  8.341e-15  8.341e-15
      // 1e-02  1.800e-13  3.625e-5  2.177e-11  2.829e-12  1.657e-14  1.135e-14
      // 1e-03  4.430e-13  2.122e-4   3.814e-8  1.831e-11  6.583e-14  6.802e-14
      // 1e-04  3.800e-11  6.053e-4   1.210e-7  1.663e-10  1.547e-13  1.117e-13
      // 1e-05   8.950e-9  4.765e-3   1.870e-7  9.665e-10  1.741e-12  4.231e-13
      // 1e-06   1.952e-7  3.108e-2   1.733e-5   1.065e-8  1.649e-11  3.719e-11
      // 1e-07   4.301e-7  2.135e-1   1.643e-5   4.137e-8  2.770e-10  1.994e-10
      // 1e-08   5.274e-8  5.634e-1   2.179e-3   2.416e-7   2.112e-9   5.444e-9
      //
      // QuadReal:
      //   eps                   32         64         96        128        160
      // 1e-01             1.461e-6  1.420e-10  2.586e-15  4.796e-17  4.354e-17
      // 1e-02             5.680e-6  2.177e-11  2.852e-12  2.578e-15  3.109e-17
      // 1e-03             1.937e-3   3.814e-8  1.827e-11  1.142e-14  2.745e-17
      // 1e-04             1.587e-2   1.210e-7  1.663e-10  1.062e-13  2.660e-17
      // 1e-05             3.396e-2   1.870e-7  9.667e-10  1.268e-12  7.498e-17
      // 1e-06             1.061e-1   1.733e-5   1.065e-8  1.061e-11  5.401e-16
      // 1e-07             3.430e-2   1.643e-5   4.138e-8  8.232e-11  9.117e-16
      // 1e-08             1.519e-2   2.179e-3   2.417e-7   1.439e-9  9.112e-15



    }

    //V[0] = test_mobility<RefReal>(R, eps_, ICIPType::Adaptive, tol*1e-14, -1, 0);
    //V[1] = test_mobility<RefReal>(R, eps_, ICIPType::Compress, tol*1e-14, tol*1e-14, 1000);
    V[2] = test_mobility<RefReal>(R, eps_, ICIPType::Precond , tol*1e-14, tol*1e-14, 100);
    //V[3] = real2quad(test_mobility<Real>(R, eps, ICIPType::Adaptive, tol, tol, 5000));
    //V[4] = real2quad(test_mobility<Real>(R, eps, ICIPType::Compress, tol, tol, 1000));
    V[5] = real2quad(test_mobility<Real>(R, eps, ICIPType::Precond , tol, tol, 100));

    //Profile::Enable(true);
    //sctl::Profile::Tic("Adap");
    //for (long i = 0; i < 0; i++) {
    //  sctl::Profile::Tic("Adap");
    //  Profile::Enable(false);
    //  V[3] = real2quad(test_mobility<Real>(R, eps, ICIPType::Adaptive, tol, tol, 5000));
    //  Profile::Enable(true);
    //  sctl::Profile::Toc();
    //}
    //sctl::Profile::Toc();
    //sctl::Profile::Tic("Precond");
    //for (long i = 0; i < 1; i++) {
    //  sctl::Profile::Tic("Precond");
    //  Profile::Enable(false);
    //  V[5] = real2quad(test_mobility<Real>(R, eps, ICIPType::Precond , tol, tol, 100));
    //  Profile::Enable(true);
    //  sctl::Profile::Toc();
    //}
    //sctl::Profile::Toc();

    Profile::print();

    for (const auto& v : V) {
      for (const auto& x : v) {
        std::cout<<std::setprecision(16)<<x<<' ';
      }
      std::cout<<'\n';
    }
    std::cout<<'\n';

    RefReal max_val = 0;
    std::cout<<"Error = ";
    for (const auto& v : V) for (const auto& x : v) max_val = std::max(max_val, fabs(x));
    for (Integer i = 0; i < 1; i++) {
      for (Integer j = 0; j < 6; j++) {
        RefReal err = 0;
        if (V[i].Dim() == V[j].Dim()) {
          for (const auto x : V[i]-V[j]) err = std::max(err, fabs(x));
        }
        std::cout<<std::setw(10)<<std::setprecision(4)<<err/max_val<<' ';
      }
      std::cout<<'\n';
    }
    std::cout<<'\n';
  }

  //const Real tol = machine_eps<Real>()*64;

  //const Real dt = 0.1;
  //const Real R = 0.75;
  //const Long Ndisc = 2;

  //Vector<Real> X(Ndisc*COORD_DIM), F(Ndisc*3);
  //{ // Set F, X
  //  const Real eps = 1e-8;

  //  X = 0;
  //  X[0] = -(R+eps/2); X[1] = 0;
  //  X[2] =  (R+eps/2); X[3] = 0;
  //  //X[0] = -(R+eps/2)*sqrt<Real>(0.5); X[1] = -(R+eps/2)*sqrt<Real>(0.5);
  //  //X[2] =  (R+eps/2)*sqrt<Real>(0.5); X[3] =  (R+eps/2)*sqrt<Real>(0.5);
  //  //X[4] =  0.0; X[5] =-0.5;

  //  F = 0;
  //  F[0] = 0; F[1] = 1; F[2] = 0;
  //  F[3] = 0; F[4] = 0; F[5] = 0;
  //  //F[6] = 0; F[7] = 0; F[8] = 0;
  //}

  //const auto real2quad = [](const Vector<Real>& v) {
  //  Vector<RefReal> w(v.Dim());
  //  for (Long i = 0; i < v.Dim(); i++) w[i] = (RefReal)v[i];
  //  return w;
  //};
  //const auto quad2real = [](const Vector<RefReal>& v) {
  //  Vector<Real> w(v.Dim());
  //  for (Long i = 0; i < v.Dim(); i++) w[i] = (Real)v[i];
  //  return w;
  //};

  //for (Long i = 0; i < 1000; i++) {
  //  const auto V0 = quad2real(DiscMobilitySolve<ElemOrder+8,RefReal>(real2quad(F), real2quad(X), R, ICIPType::Adaptive, tol*1e-3));
  //  const auto V1 = DiscMobilitySolve<ElemOrder>(F, X, R, ICIPType::Adaptive, tol);
  //  const auto V2 = DiscMobilitySolve<ElemOrder>(F, X, R, ICIPType::Compress, tol);
  //  const auto V3 = DiscMobilitySolve<ElemOrder>(F, X, R, ICIPType::Precond , tol);

  //  Real err_adap = 0, err_comp = 0, err_prec = 0, max_val = 0;
  //  for (const auto x : V0) max_val = std::max<Real>(max_val, fabs(x));
  //  for (const auto x : V0-V1) err_adap = std::max<Real>(err_adap, fabs(x));
  //  for (const auto x : V0-V2) err_comp = std::max<Real>(err_comp, fabs(x));
  //  for (const auto x : V0-V3) err_prec = std::max<Real>(err_prec, fabs(x));
  //  std::cout<<"Error (adaptive) = "<<err_adap/max_val<<'\n';
  //  std::cout<<"Error (compress) = "<<err_comp/max_val<<'\n';
  //  std::cout<<"Error (precond)  = "<<err_prec/max_val<<'\n';
  //  std::cout<<V0;

  //  for (Long j = 0; j < Ndisc; j++) {
  //    X[j*COORD_DIM+0] += V1[j*3+0] * dt;
  //    X[j*COORD_DIM+1] += V1[j*3+1] * dt;
  //  }
  //}

  Comm::MPI_Finalize();
  return 0;
}

