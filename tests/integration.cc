#include "honeycomb2/g2_weights.hpp"
#include <honeycomb2/honeycomb2.hpp>

#define Pi M_PI
#define Sqrt sqrt
#define Power pow
#define Cos cos
#define Sin sin
#define Log log

double T_test(double x1, double x2, double x3)
{
   return (1 - x1 * x1) * (1 - x2 * x2) * (1 - x3 * x3);
}

double T_test_der(double x1, double x2, double x3)
{
   return -2 * (1 - Power(x1, 2)) * (1 - Power(x2, 2)) * x3
        + 2 * (1 - Power(x1, 2)) * x2 * (1 - Power(x3, 2));
}

double T_test_6(double x1, double x2, double x3)
{
   return (1.0 - cos(T_test(x1, x2, x3))) * 2 * cos(M_PI * x2) / (x1 * x1 + x2 * x2 + x3 * x3);
}

double T_test_6_der(double x1, double x2, double x3)
{
   return (-2 * (-2 * x2 + 2 * x3) * Cos(Pi * x2)
           * (1 - Cos((1 - Power(x1, 2)) * (1 - Power(x2, 2)) * (1 - Power(x3, 2)))))
            / Power(Power(x1, 2) + Power(x2, 2) + Power(x3, 2), 2)
        + (2 * Pi * (1 - Cos((1 - Power(x1, 2)) * (1 - Power(x2, 2)) * (1 - Power(x3, 2)))) * Sin(Pi * x2))
              / (Power(x1, 2) + Power(x2, 2) + Power(x3, 2))
        + (2
           * (-2 * (1 - Power(x1, 2)) * (1 - Power(x2, 2)) * x3
              + 2 * (1 - Power(x1, 2)) * x2 * (1 - Power(x3, 2)))
           * Cos(Pi * x2) * Sin((1 - Power(x1, 2)) * (1 - Power(x2, 2)) * (1 - Power(x3, 2))))
              / (Power(x1, 2) + Power(x2, 2) + Power(x3, 2));
}

typedef double (*model_fnc_t)(double, double, double);

using integrator = Honeycomb::Integration::GaussKronrod<Honeycomb::Integration::GK_61>;

int main()
{

   std::vector<double> xBjs = {0.01, 0.05, 0.5, 0.7, 0.9};
   // std::vector<double> xBjs;
   // for (int i = 20; i > 0; i--) {
   //    xBjs.emplace_back(1.0 - i * 0.0495);
   // }

   Honeycomb::logger(Honeycomb::Logger::INFO, "Values of xBj: ");
   for (const double &x : xBjs) {
      Honeycomb::logger(Honeycomb::Logger::INFO, std::format("---- {:.6f}", x));
   }
   model_fnc_t test     = T_test_6;
   model_fnc_t test_der = T_test_6_der;

   // Exacts
   std::FILE *fp_out = std::fopen("g2.dat", "w");
   for (const double &xBj : xBjs) {
      auto fnc_external = [&](double xi) -> double {
         auto fnc_internal_1 = [&](double eta) -> double {
            return test_der(xi, -xi - eta, eta) / eta;
         };
         auto fnc_internal_2 = [&](double eta) -> double {
            return (test_der(xi, -xi - eta, eta) - test_der(eta, -xi - eta, xi)) / (eta + xi);
         };
         return (integrator::integrate(fnc_internal_1, -1, -xBj, 1.0e-10, 1.0e-10)
                 + integrator::integrate(fnc_internal_2, -xBj, 1 - xi, 1.0e-10, 1.0e-10))
              / xi;
      };
      std::fprintf(stderr, "\r\033[2K"); // Clean line
      std::fprintf(stderr, "Computing... ");
      const double g2 = integrator::integrate(fnc_external, xBj, 1, 1.0e-10, 1.0e-10);
      std::fprintf(fp_out, "%.16e\t%+.16e\n", xBj, g2);
      std::fprintf(stderr, "Done: %lf, %le\n", xBj, g2);
   }
   std::fprintf(stderr, "\n");
   std::fclose(fp_out);

   // Interpolation
   const size_t n         = 6;
   const double rmin      = 0.01;
   Honeycomb::Grid2D grid = Honeycomb::generate_compliant_Grid2D(n, {rmin, 0.15, 0.65, 1}, {9, 9, 7});

   Honeycomb::Discretization discr(grid);
   Eigen::VectorXd _fj = discr(test);

   std::FILE *fp_out_2 = std::fopen("g2_interp.dat", "w");

   Honeycomb::timer::mark begin = Honeycomb::timer::now();

   size_t sx = xBjs.size();

   Honeycomb::ThreadPool th(10);

   for (size_t i_x = 0; i_x < sx; i_x++) {
      // if (i_x <= 2) continue;
      const double &xBj = xBjs[i_x];
      std::fprintf(stderr, "Computing... ");
      Honeycomb::G2Weights weight_results(xBj, grid);

      double res = weight_results.weights.dot(_fj);

      const double g2 = res;
      std::fprintf(fp_out_2, "%.16e\t%+.16e\n", xBj, g2);
      std::fprintf(stderr, "Done: %lf, %le\n", xBj, g2);
   }
   std::fprintf(stderr, "\n");
   std::fclose(fp_out_2);
   Honeycomb::timer::mark end = Honeycomb::timer::now();

   Honeycomb::logger(Honeycomb::Logger::INFO, std::format("Full computation took: {:.4e} (ms)",
                                                          Honeycomb::timer::elapsed_ms(end, begin)));
}
