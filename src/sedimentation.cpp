#include <bie2d.hpp>
#include <fstream>
#include <random>
using namespace sctl;

constexpr Integer ElemOrder = 24;
constexpr Integer COORD_DIM = 2;
using Real = double;

std::mt19937 rng;

Vector<Real> init_figure()
{
    Vector<Real> X{0.246269,-1.01665,0.425949,-1.33609,0.836805,-0.796383,-0.526118,-1.12274,-0.242385,-0.303796,-0.616053,0.197841,0.513863,1.03774,1.19888,-0.159532,0.756344,0.37438,0.160115,-0.134131,0.431363,0.775653,0.254079,0.658974,1.18762,0.773148,0.247955,1.35082,-0.778471,-0.773944,0.715725,1.0452,-0.701955,1.08701,0.461887,-1.03104,0.0513161,-1.45635,0.0878169,0.244029,0.836175,-0.17563,0.980329,0.317331,-0.286062,0.0570924,0.48943,-0.798164,-0.0194157,0.0605511,0.543891,0.437156,0.607455,0.637777,-0.159799,-0.11945,-1.16931,0.0104415,-1.26042,0.287914,-0.140155,1.41827,-0.555323,-0.0173472,1.17629,0.0426624,0.55686,-0.434869,-0.670061,0.653058,-0.00413062,-0.25172,-0.577904,-0.926157,-1.37877,0.953169,-0.860545,0.859801,0.0353042,-0.908492,0.575603,-0.230657,-1.03196,1.06205,0.0412923,0.885524,-0.0548593,-1.08925,0.899184,-1.2346,-0.198521,-0.907329,0.273049,0.42749,-0.770991,-1.32653,0.400671,-0.00918193,-0.960612,-0.976715,-0.856385,0.0375726,-0.614396,-0.653943,0.18624,0.0661726,-0.431101,0.303565,-0.194111,0.78393,-1.07596,0.194554,-1.08433,-0.386413,-0.85164,0.305035,-0.0274634,0.493886,-0.205624,-1.25007,-0.264699,0.973195,-1.09053,0.397182,0.314259,0.222426,-0.216642,-0.504149,1.08111,-0.323646,-0.244893,0.516564,-0.0212263,-0.452996,0.680095,-0.0175777,0.31448,1.00533,-1.00053,0.701138,-0.400442,-0.622202,-0.677408,1.2934,0.0532732,0.680907,-0.651478,0.397117,1.04362,-0.543152,0.624892,0.207073,0.829553,-0.55648,0.215426,-0.778047,0.757135,-0.361525,0.681301,-1.18611,-0.308931,-1.07648,-0.475443,0.531458,1.10502,1.00762,-0.216701,-0.706149,0.376554,-0.265062,-0.483815,0.914784,-0.820882,-0.16397,-0.199785,0.248993,-1.53221,0.152421,-0.0209787,-0.655946,0.780858,-1.01035,1.03813,0.562082,0.114553,1.08898,-0.436652,-0.24844,-0.900893,0.500938,0.943841,0.00312264,0.15356,-0.554257,-0.162785,1.1476,0.355521,-0.55029,-0.730099,-0.390353};
    X *= 7.5;
    return X;
}

Vector<Real> init_two(const Real R, const Real eps)
{
    Vector<Real> X(4);
    X[0] = -(R+eps/2); X[1] = 0;
    X[2] =  (R+eps/2); X[3] = 0;
    return X;
}

Vector<Real> init_chain(const Long Ndisc, const Real R, const Real eps)
{
    Vector<Real> X(Ndisc*COORD_DIM);
    X = 0;
    Real x = 0;
    for (Long i = 0; i < Ndisc; i++) {
        X[i*COORD_DIM] = x;
        x += 2*R+eps;
    }
    return X;
}

Vector<Real> init_lattice(const Long Ndisc, const Real R, const Real eps)
{
    Vector<Real> X(Ndisc*COORD_DIM);
    X[0] = 0; X[1] = 0;
    for (Long i = 1; i < Ndisc; i++) {
        Long layer = (Long)std::round(std::sqrt(i/3.0));
        Long first = 3*layer*(layer-1) + 1;
        Long side  = (Long)std::floor((i-first)/layer);
        Long idx  = (i-first) % layer;
        X[i*COORD_DIM+0] =  layer*std::cos((side-1)*const_pi<Real>()/3) + (idx+1)*std::cos((side+1)*const_pi<Real>()/3);
        X[i*COORD_DIM+1] = -layer*std::sin((side-1)*const_pi<Real>()/3) - (idx+1)*std::sin((side+1)*const_pi<Real>()/3);
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

Vector<Real> init_random(const Long Ndisc, const Real R, const Real eps, const Real a = -1, const Real b = 1)
{
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

int main(int argc, char** argv)
{
    Comm::MPI_Init(&argc, &argv);
    rng.seed(0);

    const Comm comm = Comm::Self();
    const Real R = 0.75;
    const Real a = -4;
    const Real b =  4;

    commandline_option_start(argc, argv, "Solve the Stokes mobility problem "
        "for randomly distributed discs in sedimentation flow\nand write the disc "
        "coordinates at each time step.\n", comm);
    const bool verbose    = strtob(commandline_option(argc, argv, "--verbose",    "false", false, true,  "Verbosity", comm));
    Long Ndisc            = strtol(commandline_option(argc, argv, "--ndisc",      "3",     false, false, "Number of discs", comm), nullptr, 10);
    const Real tol        = strtod(commandline_option(argc, argv, "--tol",        "1e-12", false, false, "Discretization tolerance", comm), nullptr);
    const Real gmres_tol  = strtod(commandline_option(argc, argv, "--gmres-tol",  "1e-12", false, false, "GMRES tolerance", comm), nullptr);
    const Long gmres_iter = strtol(commandline_option(argc, argv, "--gmres-iter", "2000",  false, false, "GMRES iterations", comm), nullptr, 10);
    const Real gravity    = strtod(commandline_option(argc, argv, "--gravity",    "-1",    false, false, "Gravity", comm), nullptr);
    const Real dt         = strtod(commandline_option(argc, argv, "--dt",         "0.1",   false, false, "Time step size", comm), nullptr);
    const Long nsteps     = strtol(commandline_option(argc, argv, "--nsteps",     "100",   false, false, "Number of time steps", comm), nullptr, 10);
    const Real eps        = strtod(commandline_option(argc, argv, "--eps",        "1e-2",  false, false, "Closeness", comm), nullptr);
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
        throw std::runtime_error("Invalid initial condition.");
    }(commandline_option(argc, argv, "--init", "random", false, false, "Initial condition", comm));

    commandline_option_end(argc, argv);

    // Initial conditions
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
    }
    Ndisc = X.Dim() / COORD_DIM;

    Vector<Real> F(Ndisc*3);
    F = 0;
    // if (gravity) {
    //     for (Long i = 0; i < Ndisc; i++) {
    //         F[i*3+1] = gravity;
    //     }
    // }
    //F[0] = 1;
    F[1] = 1;

    for (const auto& x : X) {
        std::cout << std::setprecision(16) << x << ' ';
    }
    std::cout << std::endl;

    for (const auto& f : F) {
        std::cout << std::setprecision(16) << f << ' ';
    }
    std::cout << std::endl;

    Vector<Real> V;
    DiscMobility<Real,ElemOrder> disc_mobility(comm, verbose);
    std::ofstream file("sedimentation.txt");

    for (Long i = 0; i < nsteps; i++) {

        std::cout << "i = " << i << " / " << nsteps << std::endl;

        disc_mobility.Init(X, R, tol, icip_type);
        disc_mobility.Solve(V, F, Vector<Real>(), gmres_tol, gmres_iter);

        for (Long j = 0; j < Ndisc; j++) {
            X[j*COORD_DIM+0] += V[j*3+0] * dt;
            X[j*COORD_DIM+1] += V[j*3+1] * dt;
        }
        for (const auto& x : X) {
            file << std::setprecision(16) << x << ' ';
        }
        file << std::endl;
    }

    Comm::MPI_Finalize();
    return 0;
}
