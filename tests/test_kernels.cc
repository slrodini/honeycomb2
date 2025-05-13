#include <honeycomb2/honeycomb2.hpp>
#include "honeycomb_120_25_grid_points.hpp"

double T_test(double x1, double x2, double x3)
{
   return (1 - x1 * x1) * (1 - x2 * x2) * (1 - x3 * x3);
}

int main()
{
   const double Nc = 3; // NC = 1 for tests

   const size_t n         = 8;
   const double rmin      = 0.01;
   Honeycomb::Grid2D grid = Honeycomb::generate_compliant_Grid2D(n, {rmin, 0.1, 0.4, 1}, {13, 9, 7});

   Honeycomb::Kernels ker(grid, Nc, false);

   ker.ComputeTestKernel();

   Honeycomb::Discretization discr(grid);

   Eigen::VectorXd fj = discr(T_test);

   Eigen::VectorXd fj2 = ker.TestMatrix * fj;

   std::FILE *fp_3 = std::fopen("CE_single_convolution.dat", "w");

   for (Honeycomb::RnC::Triplet &x123 : honeycomb_points) {
      Honeycomb::RnC::Pair rhophi = Honeycomb::RnC::from_x123_to_rhophi(x123);

      std::fprintf(fp_3, "%.16e\t%.16e\t%.16e\t", x123[0], x123[1], x123[2]);
      std::fprintf(fp_3, "%.16e\t%.16e\n", discr.interpolate_as_weights_v3(rhophi, fj),
                   discr.interpolate_as_weights_v3(rhophi, fj2));
   }
   std::fclose(fp_3);
}