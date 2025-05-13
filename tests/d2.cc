#include "honeycomb2/discretization.hpp"
#include <honeycomb2/honeycomb2.hpp>
#include "original_model_functions.hpp"

int main()
{

   // double fnc_elapsed           = 0;
   // Honeycomb::timer::mark begin = Honeycomb::timer::now();
   // Honeycomb::timer::mark end   = Honeycomb::timer::now();

   // // fnc_elapsed                = Honeycomb::timer::elapsed_ns(end, begin);

   // const double Nc = 3; // NC = 1 for tests

   // const size_t n         = 8;
   // const double rmin      = 0.01;
   // Honeycomb::Grid2D grid = Honeycomb::generate_compliant_Grid2D(n, {rmin, 0.1, 0.4, 1}, {13, 9, 7});

   // Honeycomb::logger(Honeycomb::Logger::INFO,
   //                   std::format("Total grid size: {:d}x{:d}", grid.size, grid.size));
   // Honeycomb::logger(Honeycomb::Logger::INFO,
   //                   std::format("Total x123 size: {:d}x{:d}", grid.c_size, grid.c_size));

   // Honeycomb::Discretization discr(grid);

   // Eigen::VectorXd discr(Tu_test);
   return 0;
}