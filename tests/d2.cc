#include "honeycomb2/discretization.hpp"
#include <honeycomb2/honeycomb2.hpp>
#include "obs_weights.hpp"
#include "original_model_functions.hpp"
#include "timer.hpp"

double T_test(double x1, double x2, double x3)
{
   return (1 - x1 * x1) * (1 - x2 * x2) * (1 - x3 * x3);
}

inline double max3(double a, double b, double c)
{
   return std::max(a, std::max(b, c));
}

double T_test_2(double x1, double x2, double x3)
{
   double r = max3(fabs(x1), fabs(x2), fabs(x3));
   return (1 - x1 * x1) * (1 - x2 * x2) * (1 - x3 * x3) * log(r);
}

double DT_test_2(double x1, double x2, double x3)
{
   double r = max3(fabs(x1), fabs(x2), fabs(x3));
   return sin(M_PI * x2) * (1 - x1 * x1) * (1 - x2 * x2) * (1 - x3 * x3) / (r);
}

int main()
{

   double fnc_elapsed           = 0;
   Honeycomb::timer::mark begin = Honeycomb::timer::now();
   Honeycomb::timer::mark end   = Honeycomb::timer::now();

   // fnc_elapsed                = Honeycomb::timer::elapsed_ns(end, begin);

   // const double Nc = 3; // NC = 1 for tests

   const size_t n         = 8;
   const double rmin      = 0.001;
   Honeycomb::Grid2D grid = Honeycomb::generate_compliant_Grid2D(n, {rmin, 0.1, 0.4, 1}, {13, 9, 7});

   Honeycomb::logger(Honeycomb::Logger::INFO,
                     std::format("Total grid size: {:d}x{:d}", grid.size, grid.size));
   Honeycomb::logger(Honeycomb::Logger::INFO,
                     std::format("Total x123 size: {:d}x{:d}", grid.c_size, grid.c_size));

   Honeycomb::Discretization discr(grid);

   // Eigen::VectorXd F_test = discr(T_test);
   // Eigen::VectorXd F_test = discr(orignal_models::Tu_test);
   Eigen::VectorXd F_test = discr(T_test_2);
   // Eigen::VectorXd F_test = discr(DT_test_2);

   begin = Honeycomb::timer::now();
   Honeycomb::D2Weights d2_weights(grid, 1.0e-10);

   end         = Honeycomb::timer::now();
   fnc_elapsed = Honeycomb::timer::elapsed_ms(end, begin);
   std::cout << "Elapsed (ms): " << fnc_elapsed << std::endl;

   begin = Honeycomb::timer::now();
   Honeycomb::D2WeightsCutted d2_weights_cutted(grid, 1.0e-10);
   end         = Honeycomb::timer::now();
   fnc_elapsed = Honeycomb::timer::elapsed_ms(end, begin);
   std::cout << "Elapsed (ms): " << fnc_elapsed << std::endl;

   begin = Honeycomb::timer::now();
   Honeycomb::D1Weights d1_weights(grid, 1.0e-10);
   end         = Honeycomb::timer::now();
   fnc_elapsed = Honeycomb::timer::elapsed_ms(end, begin);
   std::cout << "Elapsed (ms): " << fnc_elapsed << std::endl;

   std::cout << std::format("{:.10f}", d2_weights.ComputeSingleQuark(F_test)) << std::endl;
   std::cout << std::format("{:.10f}", d2_weights_cutted.ComputeSingleQuark(F_test)) << std::endl;
   std::cout << std::format("{:.10f}", d2_weights_cutted.ComputeSingleQuark_NoCorrections(F_test))
             << std::endl;

   std::cout << std::format("{:.10f}", d1_weights.ComputeSingleQuark(F_test)) << std::endl;
   std::cout << std::format("{:.10f}", d1_weights.ComputeSingleQuark_NoCorrections(F_test)) << std::endl;

   return 0;
}