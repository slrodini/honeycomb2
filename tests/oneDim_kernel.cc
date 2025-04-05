#include <honeycomb2/honeycomb2.hpp>

double test(double x)
{
   return x * (1 - x);
}

int main()
{
   using integrator = Honeycomb::Integration::GaussKronrod<Honeycomb::Integration::GK_61>;

   Honeycomb::Grid grid({{0.01, 0.15, 1.0}, {10, 10}});
   Honeycomb::Discretization1D discr(grid, test);

   Eigen::MatrixXd H = Eigen::MatrixXd::Zero(grid.c_size, grid.c_size);

   for (size_t i = 0; i < grid.c_size; i++) {
      double x = grid._coord[i];
      for (size_t iP = 0; iP < grid.size; iP++) {
         auto ker = [&](double y) -> double {
            return (grid._weights[iP](y) - grid._weights[iP](x)) / (y - x);
         };
         auto supp    = grid.get_phys_support_weight_aj(iP);
         double lower = std::max(x, supp.first);

         if (supp.second <= lower) continue;

         double res = integrator::integrate(ker, lower, supp.second);

         H(i, grid._from_iw_to_ic[iP]) += res;
         if (std::fabs(grid._weights[iP](x)) > 1.0e-12) {
            H(i, grid._from_iw_to_ic[iP]) -= grid._weights[iP](x) * log((1.0 - x) / (supp.second - x));
         }
      }
   }

   Eigen::VectorXd fj2 = H * discr._fj;

   std::FILE *fp = fopen("Monodim.dat", "w");
   for (long int i = 0; i < grid.c_size_li; i++) {
      std::fprintf(fp, "%.16e\t%.16e\n", grid._coord[i], fj2(i));
   }
   std::fclose(fp);
   return 0;
}