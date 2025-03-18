#include <bie2d.hpp>
#include <fstream>
#include <random>
using namespace sctl;

constexpr Integer ElemOrder = 24;
constexpr Integer COORD_DIM = 2;
using Real = double; //QuadReal;

std::mt19937 rng;

Vector<Real> init_figure() {
  Vector<Real> X{0.246269,-1.01665,0.425949,-1.33609,0.836805,-0.796383,-0.526118,-1.12274,-0.242385,-0.303796,-0.616053,0.197841,0.513863,1.03774,1.19888,-0.159532,0.756344,0.37438,0.160115,-0.134131,0.431363,0.775653,0.254079,0.658974,1.18762,0.773148,0.247955,1.35082,-0.778471,-0.773944,0.715725,1.0452,-0.701955,1.08701,0.461887,-1.03104,0.0513161,-1.45635,0.0878169,0.244029,0.836175,-0.17563,0.980329,0.317331,-0.286062,0.0570924,0.48943,-0.798164,-0.0194157,0.0605511,0.543891,0.437156,0.607455,0.637777,-0.159799,-0.11945,-1.16931,0.0104415,-1.26042,0.287914,-0.140155,1.41827,-0.555323,-0.0173472,1.17629,0.0426624,0.55686,-0.434869,-0.670061,0.653058,-0.00413062,-0.25172,-0.577904,-0.926157,-1.37877,0.953169,-0.860545,0.859801,0.0353042,-0.908492,0.575603,-0.230657,-1.03196,1.06205,0.0412923,0.885524,-0.0548593,-1.08925,0.899184,-1.2346,-0.198521,-0.907329,0.273049,0.42749,-0.770991,-1.32653,0.400671,-0.00918193,-0.960612,-0.976715,-0.856385,0.0375726,-0.614396,-0.653943,0.18624,0.0661726,-0.431101,0.303565,-0.194111,0.78393,-1.07596,0.194554,-1.08433,-0.386413,-0.85164,0.305035,-0.0274634,0.493886,-0.205624,-1.25007,-0.264699,0.973195,-1.09053,0.397182,0.314259,0.222426,-0.216642,-0.504149,1.08111,-0.323646,-0.244893,0.516564,-0.0212263,-0.452996,0.680095,-0.0175777,0.31448,1.00533,-1.00053,0.701138,-0.400442,-0.622202,-0.677408,1.2934,0.0532732,0.680907,-0.651478,0.397117,1.04362,-0.543152,0.624892,0.207073,0.829553,-0.55648,0.215426,-0.778047,0.757135,-0.361525,0.681301,-1.18611,-0.308931,-1.07648,-0.475443,0.531458,1.10502,1.00762,-0.216701,-0.706149,0.376554,-0.265062,-0.483815,0.914784,-0.820882,-0.16397,-0.199785,0.248993,-1.53221,0.152421,-0.0209787,-0.655946,0.780858,-1.01035,1.03813,0.562082,0.114553,1.08898,-0.436652,-0.24844,-0.900893,0.500938,0.943841,0.00312264,0.15356,-0.554257,-0.162785,1.1476,0.355521,-0.55029,-0.730099,-0.390353};
  X *= 7.5;
  return X;
}

Vector<Real> init_two(const Real R, const Real eps) {
  Vector<Real> X(4);
  X[0] = -(R+eps/2); X[1] = 0;
  X[2] =  (R+eps/2); X[3] = 0;
  return X;
}

Vector<Real> init_chain(const Long Ndisc, const Real R, const Real eps) {
  Vector<Real> X(Ndisc*COORD_DIM);
  X = 0;
  Real x = 0;
  for (Long i = 0; i < Ndisc; i++) {
    X[i*COORD_DIM+1] = x;
    x += 2*R+eps;
  }
  return X;
}

Vector<Real> init_lattice(const Long Ndisc, const Real R, const Real eps) {
  Vector<Real> X(Ndisc*COORD_DIM);
  X[0] = 0; X[1] = 0;
  for (Long i = 1; i < Ndisc; i++) {
    Long layer = (Long)std::round(std::sqrt(i/3.0));
    Long first = 3*layer*(layer-1) + 1;
    Long side  = (Long)std::floor((i-first)/layer);
    Long idx  = (i-first) % layer;
    X[i*COORD_DIM+0] =  layer*cos<Real>((side-1)*const_pi<Real>()/3) + (idx+1)*cos<Real>((side+1)*const_pi<Real>()/3);
    X[i*COORD_DIM+1] = -layer*sin<Real>((side-1)*const_pi<Real>()/3) - (idx+1)*sin<Real>((side+1)*const_pi<Real>()/3);
  }

  // Rescale
  for (Real& x : X) {
    x *= 2*(R+eps/2);
  }

  Real x = X[18*COORD_DIM+0] + 2*R + 1e-6;
  Real y = X[18*COORD_DIM+1];
  X.PushBack(x);
  X.PushBack(y);

  return X;
}

Vector<Real> init_random(const Long Ndisc, const Real R, const Real eps, const Real a = -1, const Real b = 1) {
  std::uniform_real_distribution<> dis((double)a, (double)b);
  auto rand = [&]() { return dis(rng); };

  Vector<Real> X(Ndisc*COORD_DIM);
  for (Long i = 0; i < Ndisc; i++) {
    bool accept = false;
    while (!accept) {
      accept = true;
      X[i*COORD_DIM+0] = rand();
      X[i*COORD_DIM+1] = rand();
      for (Long j = 0; j < i; j++) {
        Real dx = X[i*COORD_DIM+0] - X[j*COORD_DIM+0];
        Real dy = X[i*COORD_DIM+1] - X[j*COORD_DIM+1];
        Real d = sqrt<Real>(dx*dx + dy*dy);
        if (d < 2*R+eps) {
          accept = false;
          break;
        }
      }
    }
  }
  return X;
}

std::tuple<Vector<Real>,Vector<Real>> init_test(const Long Ndisc, const Real R, const Real eps) {
  Vector<Real> X, F;
  {
    const Long N = (Long)(sqrt<Real>(Ndisc/2)+0.5);
    for (Long i = 0; i < Ndisc/2; i++) {
      const Long j0 = i%(N*2+1);
      const Long j1 = i/(N*2+1)+1;
      X.PushBack((j0-N/2) * (R+eps/2));
      X.PushBack((j1*2+j0%2) * sqrt<Real>(3)*(R+eps/2));
      F.PushBack(0);
      F.PushBack(-1);
      F.PushBack(0);
    }

    for (Long i = Ndisc/2; i < Ndisc; i++) {
      const Long j0 = (i-Ndisc/2)%(N*2+1);
      const Long j1 = (i-Ndisc/2)/(N*2+1)+1;
      X.PushBack((j0-N/2) * (R+eps/2));
      X.PushBack(-(j1*2-j0%2) * sqrt<Real>(3)*(R+eps/2));
      F.PushBack(0);
      F.PushBack(1);
      F.PushBack(0);
    }
  }
  return std::make_tuple(X, F);
}

std::tuple<Vector<Real>,Vector<Real>> init_ball(const Long Ndisc, const Real R, const Real eps) {
  Vector<Real> X(Ndisc*2), F(Ndisc*3);
  X = 0;
  F = 0;

  Long cnt = 1, k0 = 1;
  while (cnt < Ndisc) {
    for (Long k = 0; k < k0; k++) {
      for (Long i = 0; i < 6; i++) {
        Real x = -k0 + 2*k;
        Real y = -k0 * sqrt<Real>(3);
        X[cnt*2+0] = x * cos<Real>(2*const_pi<Real>()*i/6) - y * sin<Real>(2*const_pi<Real>()*i/6);
        X[cnt*2+1] = x * sin<Real>(2*const_pi<Real>()*i/6) + y * cos<Real>(2*const_pi<Real>()*i/6);
        cnt++;
      }
    }
    k0++;
  }
  X *= (R+eps);

  srand(4);
  Real F_avg = 0;
  Vector<Real> FF;
  for (Long i = 0; i < Ndisc; i++) {
    FF.PushBack((drand48()-0.5)*2);
    F_avg += FF[i];
  }
  FF -= F_avg/Ndisc;
  for (Long i = 0; i < Ndisc; i++) F[i*3+0] = FF[i];

  return std::make_tuple(X, F);
}

int main(int argc, char** argv) {
  Comm::MPI_Init(&argc, &argv);
  rng.seed(0);

  const Comm comm = Comm::Self();
  const Real R = 0.75;
  const Real a = -4;
  const Real b =  4;

  commandline_option_start(argc, argv, "Solve the Stokes mobility problem "
      "for randomly distributed discs in sedimentation flow\nand write the disc "
      "coordinates at each time step.\n", comm);
  const bool verbose         = strtob(commandline_option(argc, argv, "--verbose",    "false", false, true,  "Verbosity", comm));
  Long Ndisc                 = strtol(commandline_option(argc, argv, "--ndisc",      "3",     false, false, "Number of discs", comm), nullptr, 10);
  const Real tol             = strtod(commandline_option(argc, argv, "--tol",        "1e-14", false, false, "Discretization tolerance", comm), nullptr);
  const Real gmres_tol       = strtod(commandline_option(argc, argv, "--gmres-tol",  "1e-12", false, false, "GMRES tolerance", comm), nullptr);
  const Long gmres_iter      = strtol(commandline_option(argc, argv, "--gmres-iter", "2000",  false, false, "GMRES iterations", comm), nullptr, 10);
  const Real gravity         = strtod(commandline_option(argc, argv, "--gravity",    "-1",    false, false, "Gravity", comm), nullptr);
  const bool ts_adap         = strtob(commandline_option(argc, argv, "--ts-adap",    "false", false, true,  "Time step adaptivity", comm));
  const Long ts_order        = strtol(commandline_option(argc, argv, "--ts-order",   "5",     false, false, "Time step order", comm), nullptr, 10);
  const Real ts_tol          = strtod(commandline_option(argc, argv, "--ts-tol",     "1e-6",  false, false, "Time step tolerance", comm), nullptr);
  const Real T_end           = strtod(commandline_option(argc, argv, "--T",          "100",   false, false, "Simulation end time", comm), nullptr);
  const Real dt0             = strtod(commandline_option(argc, argv, "--dt",         "0.1",   false, false, "Time step size", comm), nullptr);
  const Real eps             = strtod(commandline_option(argc, argv, "--eps",        "1e-2",  false, false, "Closeness", comm), nullptr);
  const ICIPType icip_type = [](std::string x) {
    std::transform(x.begin(), x.end(), x.begin(), ::tolower);
    if (x == "adaptive") return ICIPType::Adaptive;
    if (x == "compress") return ICIPType::Compress;
    if (x == "precond")  return ICIPType::Precond;
    throw std::runtime_error("Invalid ICIP type.");
  }(commandline_option(argc, argv, "--icip-type", "adaptive", false, false, "ICIP type", comm));
  const std::string init = [](std::string x) {
    std::transform(x.begin(), x.end(), x.begin(), ::tolower);
    if (x == "random")  return "random";
    if (x == "two")     return "two";
    if (x == "figure")  return "figure";
    if (x == "lattice") return "lattice";
    if (x == "chain")   return "chain";
    if (x == "test")   return "test";
    if (x == "ball")   return "ball";
    throw std::runtime_error("Invalid initial condition.");
  }(commandline_option(argc, argv, "--init", "random", false, false, "Initial condition", comm));
  const std::string vis_path =        commandline_option(argc, argv, "--vis-path",   "",      false, false, "Path to save vis", comm);
  const std::string geom_fname =      commandline_option(argc, argv, "--geom-fname", "",      false, false, "Load geometry from file", comm);
  const Long start_idx       = strtol(commandline_option(argc, argv, "--start-idx",  "0",     false, false, "Starting frame index", comm), nullptr, 10);
  const bool enable_ksprecon = strtob(commandline_option(argc, argv, "--ksprecon",   "false", false, true,  "Use Krylov preconditioner", comm));
  commandline_option_end(argc, argv);

  if (!comm.Rank()) { // Print options
    std::cout<<"CMD:";
    for (Long i = 0; i < argc; i++) std::cout<<' '<<argv[i];
    std::cout<<'\n';

    std::cout<<"np         = "<<comm.Size()<<'\n';
    std::cout<<"threads    = "<<omp_get_max_threads()<<'\n';
    std::cout<<"verbose    = "<<verbose   <<'\n';
    std::cout<<"Ndisc      = "<<Ndisc     <<'\n';
    std::cout<<"tol        = "<<tol       <<'\n';
    std::cout<<"gmres_tol  = "<<gmres_tol <<'\n';
    std::cout<<"gmres_iter = "<<gmres_iter<<'\n';
    std::cout<<"gravity    = "<<gravity   <<'\n';
    std::cout<<"ts_adap    = "<<ts_adap   <<'\n';
    std::cout<<"ts_order   = "<<ts_order  <<'\n';
    std::cout<<"ts_tol     = "<<ts_tol    <<'\n';
    std::cout<<"T_end      = "<<T_end     <<'\n';
    std::cout<<"dt         = "<<dt0       <<'\n';
    std::cout<<"eps        = "<<eps       <<'\n';
    std::cout<<"icip_type  = "<<icip_type <<'\n';
    std::cout<<"init       = "<<init      <<'\n';
    std::cout<<"vis_path   = "<<vis_path  <<'\n';
    std::cout<<"geom_fname = "<<geom_fname<<'\n';
    std::cout<<"ksprecon   = "<<enable_ksprecon<<'\n';
  }

  // Initial conditions
  Vector<Real> X0, F(Ndisc*3);
  F = 0;
  if (gravity) {
    for (Long i = 0; i < Ndisc; i++) {
      F[i*3+1] = gravity;
    }
  } else {
    //F[0] = 1;
    //F[1] = 1;
    F[2] = 1;
  }

  { // Set X0 <-- coordinates + orientation
    Vector<double> geom;
    if (!geom_fname.empty()) geom.Read<double>(geom_fname.c_str());
    if (geom.Dim()) {
      const Long N = (geom.Dim()-1)/6;
      SCTL_ASSERT(R == geom[0]);
      X0.ReInit(N*3);
      F.ReInit(N*3);
      for (Long i = 0; i < N*3; i++) X0[i] = (Real)geom[i + 1];
      for (Long i = 0; i < N*3; i++) F[i] = (Real)geom[i + N*3 + 1];
    } else {
      Vector<Real> X;
      if (init == "random") {
        X = init_random(Ndisc, R, eps, a, b);
      } else if (init == "two") {
        X = init_two(R, eps);
      } else if (init == "figure") {
        X = init_figure();
      } else if (init == "lattice") {
        X = init_lattice(Ndisc, R, eps);
      } else if (init == "chain") {
        X = init_chain(Ndisc, R, eps);
      } else if (init == "test") {
        std::tie(X,F) = init_test(Ndisc, R, eps);
      } else if (init == "ball") {
        std::tie(X,F) = init_ball(Ndisc, R, eps);
      }

      const Long N = X.Dim() / COORD_DIM;
      X0.ReInit(N*(COORD_DIM+1));
      X0 = 0;
      for (Long i = 0; i < N; i++) {
        for (Long k = 0; k < COORD_DIM; k++)
          X0[i*(COORD_DIM+1)+k] = X[i*COORD_DIM+k];
      }
      SCTL_ASSERT(X0.Dim() == Ndisc*(COORD_DIM+1));
    }
  }

  DiscMobility<Real,ElemOrder> disc_mobility(comm, verbose);

  Vector<KrylovPrecond<Real>> ksprecon(enable_ksprecon ? ts_order : 0);
  Vector<DiscPanelLst<Real,ElemOrder>> ksprecon_panel_lst(enable_ksprecon ? ts_order : 0);
  auto mobility_solve = [&disc_mobility,&R,&tol,&icip_type,&F,&gmres_tol,&gmres_iter,&comm,&ksprecon,&ksprecon_panel_lst,&ts_order](Vector<Real>* V, const Vector<Real>& X, const Integer correction_idx, const Integer substep_idx) {
    //Profile::Enable(true);
    //Profile::Tic("MobilSolve", &comm, true, 1);

    const Long N = X.Dim() / (COORD_DIM+1);
    Vector<Real> X0(N * COORD_DIM);
    for (Long i = 0; i < N; i++) {
      for (Long k = 0; k < COORD_DIM; k++) {
        X0[i*COORD_DIM+k] = X[i*(COORD_DIM+1)+k];
      }
    }
    disc_mobility.Init(X0, R, tol, icip_type);
    const auto panel_lst = disc_mobility.GetPanelList();

    if (ksprecon.Dim()) {
      if (correction_idx == 0) {
        if (substep_idx == 0) {
          ksprecon[0] = ksprecon[ts_order-1];
          ksprecon_panel_lst[0] = ksprecon_panel_lst[ts_order-1];
        } else {
          ksprecon[substep_idx] = ksprecon[substep_idx-1];
          ksprecon_panel_lst[substep_idx] = ksprecon_panel_lst[substep_idx-1];
        }
      }

      if (!panel_lst.SameRefinement(ksprecon_panel_lst[substep_idx])) ksprecon[substep_idx] = KrylovPrecond<Real>();
      ksprecon_panel_lst[substep_idx].WriteVTK("vis-ksprecon");
      disc_mobility.Solve(*V, F, Vector<Real>(), gmres_tol, gmres_iter, &ksprecon[substep_idx]);
      ksprecon_panel_lst[substep_idx] = panel_lst;
    } else {
      disc_mobility.Solve(*V, F, Vector<Real>(), gmres_tol, gmres_iter, nullptr);
    }

    //{ ////////////////////// debug write VTK
    //  static Long idx = 0;
    //  disc_mobility.GetPanelList().WriteVTK(std::string("XX_")+std::to_string(idx), Vector<Real>(), comm);
    //  idx++;
    //}

    //Profile::Toc();
    //Profile::print(&comm, {"t","f","f/s","m","alloc_m","alloc_count"});
    //Profile::reset();
  };

  std::function<void(Real, Real, const Vector<Real>&)> monitor_callback = [&start_idx,&disc_mobility,&R,&tol,&icip_type,&F,&vis_path,&comm](Real t, Real dt, const Vector<Real>& X) {
    static Long idx = start_idx;

    if (!vis_path.empty()) { // Write geometry data to file
      Vector<double> geom;
      geom.PushBack((double)R);
      for (const auto x : X) geom.PushBack((double)x);
      for (const auto f : F) geom.PushBack((double)f);
      geom.Write((vis_path+std::string("/X_")+std::to_string(idx)+".geom").c_str());
    }
    if (!vis_path.empty()) { // Write vtu file
      const Long M = 200;
      const Long N = X.Dim() / (COORD_DIM+1);
      VTUData vtu_data;
      Real Xc[2];
      { // set vtk center
        Real XX0[3] = {0,0,0};
        for (Long i = 0; i < N; i++) {
          Real XX[3] = {0,0,0};
          for (Long j = 0; j < N; j++) {
            const Real dX[2] = {X[i*(COORD_DIM+1)+0]-X[j*(COORD_DIM+1)+0], X[i*(COORD_DIM+1)+1]-X[j*(COORD_DIM+1)+1]};
            const Real R2 = dX[0]*dX[0] + dX[1]*dX[1];
            const Real Rinv = (R2>0 ? 1/sqrt<Real>(R2) : 0 );
            XX[0] += X[i*(COORD_DIM+1)+0] * Rinv;
            XX[1] += X[i*(COORD_DIM+1)+1] * Rinv;
            XX[2] += Rinv;
          }
          //if (XX[2] > XX0[2]) for (Long k = 0; k < 3; k++) XX0[k] = XX[k];
          for (Long k = 0; k < 3; k++) XX0[k] += XX[k]*pow<Real>(XX[2],16);
        }
        for (Long k = 0; k < 2; k++) Xc[k] = XX0[k]/XX0[2];
      }
      for (Long i = 0; i < N; i++) {
        const Real x = X[i*(COORD_DIM+1)+0]-Xc[0];
        const Real y = X[i*(COORD_DIM+1)+1]-Xc[1];
        const auto f = F.begin() + i*(COORD_DIM+1);
        vtu_data.connect.PushBack(vtu_data.coord.Dim()/3);
        vtu_data.coord.PushBack((float)x);
        vtu_data.coord.PushBack((float)y);
        vtu_data.coord.PushBack((float)0);
        for (Long k = 0; k < 3; k++) vtu_data.value.PushBack((float)f[k]);
        for (Long j = 0; j < M; j++) {
          const Real theta = 2*const_pi<Real>()*j/(M-1) + X[i*(COORD_DIM+1)+2]/R;
          vtu_data.connect.PushBack(vtu_data.coord.Dim()/3);
          vtu_data.coord.PushBack((float)(x + R*cos<Real>(theta)));
          vtu_data.coord.PushBack((float)(y + R*sin<Real>(theta)));
          vtu_data.coord.PushBack((float)0);
          for (Long k = 0; k < 3; k++) vtu_data.value.PushBack((float)f[k]);
        }
        vtu_data.offset.PushBack(vtu_data.connect.Dim());
        vtu_data.types.PushBack(4);
      }
      vtu_data.WriteVTK(vis_path+"/X_"+std::to_string(idx), comm);

      //const Long N = X.Dim() / (COORD_DIM+1);
      //Vector<Real> X0(N * COORD_DIM);
      //for (Long i = 0; i < N; i++) {
      //  for (Long k = 0; k < COORD_DIM; k++) {
      //    X0[i*COORD_DIM+k] = X[i*(COORD_DIM+1)+k];
      //  }
      //}
      //disc_mobility.Init(X0, R, tol, icip_type);
      //disc_mobility.GetPanelList().WriteVTK(std::string("vis/X_")+std::to_string(idx), Vector<Real>(), comm);
    }
    const Real min_d = [&X,&R](){
      Real min_d = R;
      for (Long i = 0; i < X.Dim()/3; i++) {
        for (Long j = 0; j < X.Dim()/3; j++) {
          if (i != j) {
            const Real d = sqrt<Real>( (X[i*3+0]-X[j*3+0])*(X[i*3+0]-X[j*3+0]) + (X[i*3+1]-X[j*3+1])*(X[i*3+1]-X[j*3+1]) ) - 2*R;
            min_d = std::min<Real>(min_d, d);
          }
        }
      }
      return min_d;
    }();

    std::cout<<"Frame-idx = "<<idx<<"    T = "<<t<<"    dt = "<<dt;
    printf("      d_min = %.15f\n\n\n", (double)min_d);
    idx++;
  };

  Vector<Real> X = X0;
  monitor_callback(0, dt0, X);
  if (ts_adap) { // adaptive time-stepping
    SDC<Real> time_step(ts_order, comm);
    time_step.AdaptiveSolve(&X, dt0, T_end, X0, mobility_solve, ts_tol, &monitor_callback, true);
  } else {
    if (ts_order == 1) { // non-adaptive, first-order time-stepping
      Real dt = dt0;
      Vector<Real> V(X.Dim());
      for (Real t = 0; t < T_end; t += dt) {
        const Real dt_ = std::min<Real>(dt, T_end-t);
        mobility_solve(&V, X, 0, 0);
        if (V.Dim()) {
          X += V * dt_;
          monitor_callback(t+dt_, dt_, X);
        } else {
          t -= dt;
          dt *= 0.5;
        }
      }
    } else if (ts_order > 1) { // non-adaptive, high-order time-stepping
      Real dt = dt0;
      Vector<Real> X_;
      SDC<Real> time_step(ts_order, comm);
      for (Real t = 0; t < T_end; t += dt) {
        const Real dt_ = std::min<Real>(dt, T_end-t);
        time_step(&X_, dt_, X, mobility_solve, time_step.Order()*2, ts_tol*dt_/T_end);
        if (X_.Dim()) {
          X = X_;
          monitor_callback(t+dt_, dt_, X);
        } else {
          t -= dt;
          dt *= 0.5;
        }
      }
    } else SCTL_ASSERT(ts_order>0);
  }

  Comm::MPI_Finalize();
  return 0;
}
