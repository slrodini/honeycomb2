#include <honeycomb2/honeycomb2.hpp>
#include "honeycomb_120_25_grid_points.hpp"

double T_test(double x1, double x2, double x3)
{
   return (1 - x1 * x1) * (1 - x2 * x2) * (1 - x3 * x3);
}

double T_test_Hu(double x1, double x2, double x3)
{
   double r = std::max(std::fabs(x1), std::max(std::fabs(x2), std::fabs(x3)));
   return sin(x2 * M_PI) * (1 - x1 * x1) * (1 - x2 * x2) * (1 - x3 * x3) / (sqrt(r));
}

void evolve_chiral_odd();

int main()
{

   const double Nc = 3; // NC = 1 for tests

   const size_t n         = 6;
   const double rmin      = 0.01;
   Honeycomb::Grid2D grid = Honeycomb::generate_compliant_Grid2D(n, {rmin, 0.15, 0.65, 1}, {9, 9, 7});

   // const size_t n         = 13;
   // const double rmin      = 0.01;
   // Honeycomb::Grid2D grid = Honeycomb::generate_compliant_Grid2D(n, {rmin, 0.15, 0.65, 1}, {13, 13, 11});

   double fnc_elapsed = 0;
   Honeycomb::logger(Honeycomb::Logger::INFO, std::format("Total grid size: {:d}x{:d}", grid.size, grid.size));

   Honeycomb::timer::mark begin = Honeycomb::timer::now();
   Honeycomb::Kernels kers(grid, Nc);
   Honeycomb::timer::mark end  = Honeycomb::timer::now();
   fnc_elapsed                += Honeycomb::timer::elapsed_ns(end, begin);

   Honeycomb::logger(Honeycomb::Logger::INFO, std::format("Discretization time: {:+.6e} ms.", fnc_elapsed * 1.0e-6));

   for (long int i = 0; i < kers.H_NS.rows(); i++) {
      for (long int j = 0; j < kers.H_NS.cols(); j++) {
         if (std::isnan(kers.H_NS(i, j)) || std::isinf(kers.H_NS(i, j)))
            Honeycomb::logger(Honeycomb::Logger::WARNING, std::format("NAN or INF in kernels! {:d}, {:d}", i, j));
      }
   }

   Honeycomb::Discretization discr(grid);
   Eigen::VectorXd fj    = discr(T_test);
   Eigen::VectorXd fj_Hu = discr(T_test_Hu);

   Eigen::VectorXd fj2 = kers.H_NS * fj;

   std::FILE *fp = std::fopen("ConvCheck.dat", "w");
   for (long int i = 0; i < grid.c_size_li; i++) {
      std::fprintf(fp, "%.16e\t%.16e\t%.16e\t%.16e\t%.16e\n", grid._x123[i].v[0], grid._x123[i].v[1],
                   grid._x123[i].v[2], fj(i), fj2(i));
   }
   std::fclose(fp);

   return 0;
}

void evolve_chiral_odd()
{
   const double Nc = 3;

   const size_t n         = 6;
   const double rmin      = 0.01;
   Honeycomb::Grid2D grid = Honeycomb::generate_compliant_Grid2D(n, {rmin, 0.15, 0.65, 1}, {9, 9, 7});

   // const size_t n         = 13;
   // const double rmin      = 0.01;
   // Honeycomb::Grid2D grid = Honeycomb::generate_compliant_Grid2D(n, {rmin, 0.15, 0.65, 1}, {13, 13, 11});

   double fnc_elapsed = 0;
   Honeycomb::logger(Honeycomb::Logger::INFO, std::format("Total grid size: {:d}x{:d}", grid.size, grid.size));

   Honeycomb::timer::mark begin  = Honeycomb::timer::now();
   Eigen::MatrixXd H_CO          = Honeycomb::get_CO_kernel(grid, Nc);
   Honeycomb::timer::mark end    = Honeycomb::timer::now();
   fnc_elapsed                  += Honeycomb::timer::elapsed_ns(end, begin);

   Honeycomb::logger(Honeycomb::Logger::INFO, std::format("Discretization time: {:+.6e} ms.", fnc_elapsed * 1.0e-6));

   Honeycomb::Discretization discr(grid);
   Eigen::VectorXd fj    = discr(T_test);
   Eigen::VectorXd fj_Hu = discr(T_test_Hu);

   const double Q02 = 1.0;
   const double Qf2 = 10000.0;

   const double t0 = log(Q02);
   const double tf = log(Qf2);

   const double dt_2 = 0.01; // ! Unused !

   const double pref = -1.0; // Evolution Equations are df/dt = - as Hxf

   runge_kutta::GenericRungeKutta<Eigen::MatrixXd, Eigen::VectorXd, 13> evolver_CO(
       H_CO, fj, runge_kutta::DOPRI8,
       [](double t) {
          return 1.0 / (11.0 * (t + 3.0));
       },
       -1.0, t0, dt_2,
       [](double t, Eigen::MatrixXd &K, Eigen::VectorXd &S) {
          (void)t;
          (void)K;
          (void)S;
          return;
       },
       fj.size());

   fnc_elapsed = 0;
   begin       = Honeycomb::timer::now();
   evolver_CO({log(10000)}, 40);
   end          = Honeycomb::timer::now();
   fnc_elapsed += Honeycomb::timer::elapsed_ns(end, begin);

   auto fj_fin = evolver_CO.GetSolution();

   Honeycomb::logger(Honeycomb::Logger::INFO, std::format("Evolution time: {:+.6e} ms.", fnc_elapsed * 1.0e-6));

   fnc_elapsed = 0;
   begin       = Honeycomb::timer::now();
   evolver_CO.reset(fj_Hu, t0);
   evolver_CO({log(10000)}, 40);
   end          = Honeycomb::timer::now();
   fnc_elapsed += Honeycomb::timer::elapsed_ns(end, begin);

   Honeycomb::logger(Honeycomb::Logger::INFO, std::format("Evolution time: {:+.6e} ms.", fnc_elapsed * 1.0e-6));

   auto fj_Hu_fin = evolver_CO.GetSolution();

   std::FILE *fp_2 = std::fopen("CO_evo_check.dat", "w");

   for (Honeycomb::RnC::Triplet &x123 : honeycomb_points) {
      Honeycomb::RnC::Pair rhophi = Honeycomb::RnC::from_x123_to_rhophi(x123);
      std::fprintf(fp_2, "%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\n", x123.v[0], x123.v[1], x123.v[2],
                   discr.interpolate_as_weights_v3(rhophi, fj), discr.interpolate_as_weights_v3(rhophi, fj_fin),
                   discr.interpolate_as_weights_v3(rhophi, fj_Hu), discr.interpolate_as_weights_v3(rhophi, fj_Hu_fin));
   }
   std::fclose(fp_2);
}