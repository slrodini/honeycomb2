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

void test_save_and_laod();
void test_convolution();
void evolve_chiral_odd();

int main()
{

   // evolve_chiral_odd();
   // test_convolution();

   test_save_and_laod();

   return 0;
}

void test_save_and_laod()
{
   const double Nc = 3;

   const size_t n         = 6;
   const double rmin      = 0.01;
   Honeycomb::Grid2D grid = Honeycomb::generate_compliant_Grid2D(n, {rmin, 0.15, 0.65, 1}, {9, 9, 7});

   double fnc_elapsed = 0;
   Honeycomb::logger(Honeycomb::Logger::INFO, std::format("Total grid size: {:d}x{:d}", grid.size, grid.size));
   Honeycomb::logger(Honeycomb::Logger::INFO, std::format("Total x123 size: {:d}x{:d}", grid.c_size, grid.c_size));

   // Discretize the kernels
   Honeycomb::timer::mark begin = Honeycomb::timer::now();
   Honeycomb::Kernels kers(grid, Nc);
   Honeycomb::timer::mark end = Honeycomb::timer::now();
   fnc_elapsed                = Honeycomb::timer::elapsed_ns(end, begin);
   Honeycomb::logger(Honeycomb::Logger::INFO, std::format("Discretization time: {:+.6e} ms.", fnc_elapsed * 1.0e-6));

   // Save the discretized kernels
   begin = Honeycomb::timer::now();
   Honeycomb::save_kernels<cereal::PortableBinaryOutputArchive>(kers, "check.cereal");
   end         = Honeycomb::timer::now();
   fnc_elapsed = Honeycomb::timer::elapsed_ns(end, begin);
   Honeycomb::logger(Honeycomb::Logger::INFO, std::format("Saving time: {:+.6e} ms.", fnc_elapsed * 1.0e-6));

   // Load the kernels
   begin = Honeycomb::timer::now();
   Honeycomb::Kernels kers_loaded =
       Honeycomb::load_kernels<cereal::PortableBinaryInputArchive>("check.cereal", grid, Nc);
   end         = Honeycomb::timer::now();
   fnc_elapsed = Honeycomb::timer::elapsed_ns(end, begin);
   Honeycomb::logger(Honeycomb::Logger::INFO, std::format("Loading time: {:+.6e} ms.", fnc_elapsed * 1.0e-6));

   for (long int i = 0; i < kers.H_NS.rows(); i++) {
      for (long int j = 0; j < kers.H_NS.cols(); j++) {
         if (std::fabs(kers.H_NS(i, j) - kers_loaded.H_NS(i, j)) > 2.0e-16) {
            Honeycomb::logger(Honeycomb::Logger::ERROR, "H_NS is not save/loaded correctly.");
         }
         if (std::fabs(kers.H_d13(i, j) - kers_loaded.H_d13(i, j)) > 2.0e-16) {
            Honeycomb::logger(Honeycomb::Logger::ERROR, "H_d13 is not save/loaded correctly.");
         }
         if (std::fabs(kers.H_gg_p(i, j) - kers_loaded.H_gg_p(i, j)) > 2.0e-16) {
            Honeycomb::logger(Honeycomb::Logger::ERROR, "H_gg_p is not save/loaded correctly.");
         }
         if (std::fabs(kers.H_gg_m(i, j) - kers_loaded.H_gg_m(i, j)) > 2.0e-16) {
            Honeycomb::logger(Honeycomb::Logger::ERROR, "H_gg_m is not save/loaded correctly.");
         }
         if (std::fabs(kers.H_qg_p(i, j) - kers_loaded.H_qg_p(i, j)) > 2.0e-16) {
            Honeycomb::logger(Honeycomb::Logger::ERROR, "H_qg_p is not save/loaded correctly.");
         }
         if (std::fabs(kers.H_qg_m(i, j) - kers_loaded.H_qg_m(i, j)) > 2.0e-16) {
            Honeycomb::logger(Honeycomb::Logger::ERROR, "H_qg_m is not save/loaded correctly.");
         }
         if (std::fabs(kers.H_gq_p(i, j) - kers_loaded.H_gq_p(i, j)) > 2.0e-16) {
            Honeycomb::logger(Honeycomb::Logger::ERROR, "H_gq_p is not save/loaded correctly.");
         }
         if (std::fabs(kers.H_gq_m(i, j) - kers_loaded.H_gq_m(i, j)) > 2.0e-16) {
            Honeycomb::logger(Honeycomb::Logger::ERROR, "H_gq_m is not save/loaded correctly.");
         }
      }
   }
}

void test_convolution()
{
   const double Nc = 3; // NC = 1 for tests

   // const size_t n         = 6;
   // const double rmin      = 0.01;
   // Honeycomb::Grid2D grid = Honeycomb::generate_compliant_Grid2D(n, {rmin, 0.15, 0.65, 1}, {9, 9, 7});

   const size_t n         = 13;
   const double rmin      = 0.01;
   Honeycomb::Grid2D grid = Honeycomb::generate_compliant_Grid2D(n, {rmin, 0.15, 0.65, 1}, {13, 13, 11});

   double fnc_elapsed = 0;
   Honeycomb::logger(Honeycomb::Logger::INFO, std::format("Total grid size: {:d}x{:d}", grid.size, grid.size));
   Honeycomb::logger(Honeycomb::Logger::INFO, std::format("Total x123 size: {:d}x{:d}", grid.c_size, grid.c_size));

   Honeycomb::timer::mark begin = Honeycomb::timer::now();
   Honeycomb::Kernels kers(grid, Nc);
   Honeycomb::timer::mark end  = Honeycomb::timer::now();
   fnc_elapsed                += Honeycomb::timer::elapsed_ns(end, begin);

   Honeycomb::logger(Honeycomb::Logger::INFO, std::format("Discretization time: {:+.6e} ms.", fnc_elapsed * 1.0e-6));

   // const auto &KK = kers.H_NS;
   const auto &KK = kers.H_test_1;

   for (long int i = 0; i < KK.rows(); i++) {
      for (long int j = 0; j < KK.cols(); j++) {
         if (std::isnan(KK(i, j)) || std::isinf(KK(i, j)))
            Honeycomb::logger(Honeycomb::Logger::WARNING, std::format("NAN or INF in kernels! {:d}, {:d}", i, j));
      }
   }

   const auto &KK2 = kers.H_test_2;
   // Eigen::MatrixXd Diff = KK2 - KK;

   // for (long int i = 0; i < KK.rows(); i++) {
   //    for (long int j = 0; j < KK.cols(); j++) {
   //       if (std::fabs(Diff(i, j)) > 1.0e-6)
   //          std::cout << std::format("{:d}\t{:d}\t{:.10e}\t{:.10e}\t{:.10e}\t{:.10e}\t{:.10e}", i, j,
   //                                   grid._x123[i].v[0], grid._x123[i].v[1], grid._x123[i].v[2], KK(i, j), KK2(i, j))
   //                    << std::endl;
   //    }
   // }

   Honeycomb::Discretization discr(grid);
   Eigen::VectorXd fj    = discr(T_test);
   Eigen::VectorXd fj_Hu = discr(T_test_Hu);

   Eigen::VectorXd fj2 = KK * fj;
   Eigen::VectorXd fj3 = KK2 * fj;

   std::FILE *fp = std::fopen("ConvCheck.dat", "w");
   for (long int i = 0; i < grid.c_size_li; i++) {
      std::fprintf(fp, "%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\n", grid._x123[i].v[0], grid._x123[i].v[1],
                   grid._x123[i].v[2], fj(i), fj2(i), fj3(i));
   }
   std::fclose(fp);

   double m1 = 0, m2 = 0;

   for (long int i = 0; i < grid.c_size_li; i++) {
      // x132
      auto rhophi = Honeycomb::RnC::from_x123_to_rhophi(
          Honeycomb::RnC::Triplet(-grid._x123[i].v[0], -grid._x123[i].v[1], -grid._x123[i].v[2]));

      const double diff_1 = fj(i) - discr.interpolate_as_weights_v3(rhophi, fj);
      const double diff_2 = fj2(i) + discr.interpolate_as_weights_v3(rhophi, fj2);

      m1 = std::max(m1, std::fabs(diff_1));
      m2 = std::max(m2, std::fabs(diff_2));
   }
   std::fprintf(stderr, "%.16e\t%.16e\n", m1, m2);
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