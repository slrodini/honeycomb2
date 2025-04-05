#include <honeycomb2/honeycomb2.hpp>

double T_test(double x1, double x2, double x3)
{
   return (1 - x1 * x1) * (1 - x2 * x2) * (1 - x3 * x3);
}

int main()
{

   // const size_t n    = 6;
   // const double rmin = 0.01;

   // Honeycomb::Grid2D grid = Honeycomb::generate_compliant_Grid2D(n, {rmin, 0.15, 0.65, 1}, {9, 9, 7});

   const size_t n    = 13;
   const double rmin = 0.01;

   Honeycomb::Grid2D grid = Honeycomb::generate_compliant_Grid2D(n, {rmin, 0.15, 0.65, 1}, {13, 13, 11});

   double fnc_elapsed = 0;
   Honeycomb::logger(Honeycomb::Logger::INFO, std::format("Total grid size: {:d}x{:d}", grid.size, grid.size));

   Honeycomb::timer::mark begin = Honeycomb::timer::now();
   Honeycomb::Kernels kers(grid, 3);
   Honeycomb::timer::mark end  = Honeycomb::timer::now();
   fnc_elapsed                += Honeycomb::timer::elapsed_ns(end, begin);

   Honeycomb::logger(Honeycomb::Logger::INFO, std::format("Discretization time: {:+.6e} ms.", fnc_elapsed * 1.0e-6));

   for (long int i = 0; i < kers.H_NS.rows(); i++) {
      for (long int j = 0; j < kers.H_NS.cols(); j++) {
         if (std::isnan(kers.H_NS(i, j)) || std::isinf(kers.H_NS(i, j)))
            Honeycomb::logger(Honeycomb::Logger::WARNING, std::format("NAN or INF in kernels! {:d}, {:d}", i, j));
      }
   }

   double max = 0;
   for (long int i = 0; i < kers.H_NS.rows(); i++) {
      for (long int j = 0; j < kers.H_NS.cols(); j++) {
         max = std::max(max, std::fabs(kers.H_plus_12_v1(i, j) - kers.H_plus_12_v2(i, j)));
      }
   }
   Honeycomb::logger(Honeycomb::Logger::INFO, std::format("Hplus12 max diff: {:.16e}", max));

   Honeycomb::Discretization discr(grid);
   Eigen::VectorXd fj = discr(T_test);

   Eigen::VectorXd fj2 = kers.H_NS * fj;

   std::FILE *fp = std::fopen("ConvCheck.dat", "w");
   for (long int i = 0; i < grid.c_size_li; i++) {
      std::fprintf(fp, "%.16e\t%.16e\t%.16e\t%.16e\t%.16e\n", grid._x123[i].v[0], grid._x123[i].v[1],
                   grid._x123[i].v[2], fj(i), fj2(i));
   }
   std::fclose(fp);

   return 0;
}
