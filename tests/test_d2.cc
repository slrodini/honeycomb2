#include <honeycomb2/honeycomb2.hpp>

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

   const size_t n    = 8;
   const double rmin = 0.0001;
   Honeycomb::Grid2D grid
       = Honeycomb::generate_compliant_Grid2D(n, {rmin, 0.01, 0.1, 0.4, 1}, {13, 9, 9, 7});

   Honeycomb::logger(Honeycomb::Logger::INFO,
                     std::format("Total grid size: {:d}x{:d}", grid.size, grid.size));
   Honeycomb::logger(Honeycomb::Logger::INFO,
                     std::format("Total x123 size: {:d}x{:d}", grid.c_size, grid.c_size));

   Honeycomb::Discretization discr(grid);

   Eigen::VectorXd F_test_1 = discr(T_test);
   Eigen::VectorXd F_test_2 = discr(T_test_2);
   Eigen::VectorXd F_test_3 = discr(DT_test_2);

   begin = Honeycomb::timer::now();
   Honeycomb::D2Weights d2_weights(grid, true, 1.0e-10);

   end         = Honeycomb::timer::now();
   fnc_elapsed = Honeycomb::timer::elapsed_ms(end, begin);
   std::cout << "Elapsed (ms): " << fnc_elapsed << std::endl;

   begin = Honeycomb::timer::now();
   Honeycomb::D2WeightsCutted d2_weights_cutted(grid, true, 1.0e-10);
   end         = Honeycomb::timer::now();
   fnc_elapsed = Honeycomb::timer::elapsed_ms(end, begin);
   std::cout << "Elapsed (ms): " << fnc_elapsed << std::endl;

   begin = Honeycomb::timer::now();
   Honeycomb::ELTWeights elt_weights(grid, true, 1.0e-10);
   end         = Honeycomb::timer::now();
   fnc_elapsed = Honeycomb::timer::elapsed_ms(end, begin);
   std::cout << "Elapsed (ms): " << fnc_elapsed << std::endl;

   double res_d2_1  = 47.0 / 160.0;
   double res_elt_1 = 913.0 / 3150.0;

   double res_d2_2  = -(949.0 / 3840.0);
   double res_elt_2 = -(69634.0 / 165375.0);

   if (!Honeycomb::is_near(res_d2_1, d2_weights.ComputeSingleQuark(F_test_1))) {
      Honeycomb::logger(Honeycomb::Logger::ERROR, "Failed d2 for T_test");
   }
   if (!Honeycomb::is_near(res_elt_1, elt_weights.ComputeSingleQuark(F_test_1), 1.0e-9)) {
      Honeycomb::logger(Honeycomb::Logger::ERROR, "Failed elt for T_test");
   }

   if (!Honeycomb::is_near(res_d2_2, d2_weights.ComputeSingleQuark(F_test_2), 1.0e-5)) {
      Honeycomb::logger(Honeycomb::Logger::ERROR, "Failed d2 for T_test_2");
   }
   if (!Honeycomb::is_near(res_elt_2, elt_weights.ComputeSingleQuark(F_test_2), 1.0e-4)) {
      Honeycomb::logger(Honeycomb::Logger::ERROR, "Failed elt for T_test_2");
   }

   if (!Honeycomb::is_near(0.0, d2_weights.ComputeSingleQuark(F_test_3), 1.0e-10)) {
      Honeycomb::logger(Honeycomb::Logger::ERROR, "Failed d2 for DT_test");
   }
   if (!Honeycomb::is_near(0.7033003506114330, elt_weights.ComputeSingleQuark(F_test_3), 1.0e-10)) {
      Honeycomb::logger(Honeycomb::Logger::ERROR, "Failed elt for DT_test");
   }

   return 0;
}