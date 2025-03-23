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
   return -2 * (1 - Power(x1, 2)) * (1 - Power(x2, 2)) * x3 + 2 * (1 - Power(x1, 2)) * x2 * (1 - Power(x3, 2));
}

double T_test_6(double x1, double x2, double x3)
{
   return (1.0 - cos(T_test(x1, x2, x3))) * 2 * cos(M_PI * x2) / (x1 * x1 + x2 * x2 + x3 * x3);
}

double T_test_6_der(double x1, double x2, double x3)
{
   return (-2 * (-2 * x2 + 2 * x3) * Cos(Pi * x2) *
           (1 - Cos((1 - Power(x1, 2)) * (1 - Power(x2, 2)) * (1 - Power(x3, 2))))) /
              Power(Power(x1, 2) + Power(x2, 2) + Power(x3, 2), 2) +
          (2 * Pi * (1 - Cos((1 - Power(x1, 2)) * (1 - Power(x2, 2)) * (1 - Power(x3, 2)))) * Sin(Pi * x2)) /
              (Power(x1, 2) + Power(x2, 2) + Power(x3, 2)) +
          (2 * (-2 * (1 - Power(x1, 2)) * (1 - Power(x2, 2)) * x3 + 2 * (1 - Power(x1, 2)) * x2 * (1 - Power(x3, 2))) *
           Cos(Pi * x2) * Sin((1 - Power(x1, 2)) * (1 - Power(x2, 2)) * (1 - Power(x3, 2)))) /
              (Power(x1, 2) + Power(x2, 2) + Power(x3, 2));
}

typedef double (*model_fnc_t)(double, double, double);

int main()
{
   using integrator = Honeycomb::Integration::GaussKronrod<Honeycomb::Integration::GK_61>;

   std::vector<double> xBjs = {0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 0.9};

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
         return (integrator::integrate(fnc_internal_1, -1, -xBj, 1.0e-5, 1.0e-6).first +
                 integrator::integrate(fnc_internal_2, -xBj, 1 - xi, 1.0e-5, 1.0e-6).first) /
                xi;

         // return (integrator::integrate(fnc_internal_1, -1, -xBj, 1.0e-5, 1.0e-6).first) / xi;
      };
      std::fprintf(stderr, "\r\033[2K"); // Clean line
      std::fprintf(stderr, "Computing... ");
      const double g2 = integrator::integrate(fnc_external, xBj, 1, 1.0e-5, 1.0e-6).first;
      std::fprintf(fp_out, "%.16e\t%+.16e\n", xBj, g2);
      std::fprintf(stderr, "Done: %lf, %le", xBj, g2);
   }
   std::fprintf(stderr, "\n");
   std::fclose(fp_out);

   // Interpolation
   const size_t n    = 12;
   const double rmin = 0.01;

   Honeycomb::Grid2D grid = Honeycomb::generate_compliant_Grid2D(n, {rmin, 0.33, 0.66, 1}, {12, 12, 12});
   Honeycomb::Discretization f(grid, test);

   fp_out = std::fopen("g2_interp.dat", "w");

   std::vector<std::map<size_t, double>> x_weight_results;

   Honeycomb::timer::mark begin = Honeycomb::timer::now();
   for (const double &xBj : xBjs) {

      std::fprintf(stderr, "Computing... \n");
      std::map<size_t, double> weight_results;

      double res = 0;
      //  Sector 0, only first integral
      for (size_t k = grid.grid_angle._delim_indexes[0]; k < grid.grid_angle._delim_indexes[1]; k++) {
         for (size_t j = 0; j < grid.grid_radius.size; j++) {
            auto r_supp          = grid.grid_radius.get_support_weight_aj(j);
            const double eta_min = std::max(grid.grid_radius._d_info.to_phys_space(r_supp.first), xBj);
            const double eta_max = std::min(grid.grid_radius._d_info.to_phys_space(r_supp.second), 1.0);
            auto index           = f.get_flatten_index(j, k);

            auto fnc_external = [&](double eta_abs) -> double {
               const double eta    = -eta_abs;
               auto fnc_internal_1 = [&](double xi) -> double {
                  auto rhophi = Honeycomb::RnC::from_x123_to_rhophi(xi, -eta - xi, eta);
                  return f._dw_dx3[index](rhophi) / xi;
               };

               return integrator::integrate(fnc_internal_1, xBj, -eta, 1.0e-5, 1.0e-6).first / eta;
            };
            const double tmp = integrator::integrate(fnc_external, eta_min, eta_max, 1.0e-5, 1.0e-6).first;
            // Works both if key does and does not exist.
            weight_results[index] += tmp;
            // res += tmp * f._fj(index);
            std::fprintf(stderr, "\r\033[2K"); // Clean line
            std::fprintf(stderr, "Done: %zu, %zu", j, k);
         }
      }

      // Sector 3, only third integral (the one with switched arguments)
      for (size_t k = grid.grid_angle._delim_indexes[3]; k < grid.grid_angle._delim_indexes[4]; k++) {
         for (size_t j = 0; j < grid.grid_radius.size; j++) {
            auto r_supp         = grid.grid_radius.get_support_weight_aj(j);
            const double xi_min = std::max(grid.grid_radius._d_info.to_phys_space(r_supp.first), xBj);
            const double xi_max = std::min(grid.grid_radius._d_info.to_phys_space(r_supp.second), 1.0);
            auto index          = f.get_flatten_index(j, k);

            auto fnc_external = [&](double xi) -> double {
               auto fnc_internal_1 = [&](double eta) -> double {
                  // Switched arguments, and mind the global `-` sign
                  auto rhophi = Honeycomb::RnC::from_x123_to_rhophi(eta, -eta - xi, xi);
                  return -f._dw_dx3[index](rhophi) / (eta + xi);
               };

               return (integrator::integrate(fnc_internal_1, -xBj, 0, 1.0e-5, 1.0e-6).first) / xi;
            };
            const double tmp = integrator::integrate(fnc_external, xi_min, xi_max, 1.0e-5, 1.0e-6).first;
            // Works both if key does and does not exist.
            weight_results[index] += tmp;
            // res += tmp * f._fj(index);
            std::fprintf(stderr, "\r\033[2K"); // Clean line
            std::fprintf(stderr, "Done: %zu, %zu", j, k);
         }
      }

      // Sector 4, second and third integrals, written as integrals over t
      for (size_t k = grid.grid_angle._delim_indexes[4]; k < grid.grid_angle._delim_indexes[5]; k++) {
         for (size_t j = 0; j < grid.grid_radius.size; j++) {
            auto r_supp        = grid.grid_radius.get_support_weight_aj(j);
            const double t_min = std::max(grid.grid_radius._d_info.to_phys_space(r_supp.first), xBj);
            const double t_max = std::min(grid.grid_radius._d_info.to_phys_space(r_supp.second), 1.0);
            auto index         = f.get_flatten_index(j, k);

            auto fnc_external = [&](double t_abs) -> double {
               const double t      = -t_abs;
               auto fnc_internal_1 = [&](double xi) -> double {
                  auto rhophi_1 = Honeycomb::RnC::from_x123_to_rhophi(xi, t, -xi - t);
                  auto rhophi_2 = Honeycomb::RnC::from_x123_to_rhophi(-xi - t, t, xi);
                  return (f._dw_dx3[index](rhophi_1) - f._dw_dx3[index](rhophi_2)) / xi;
               };

               // yes, |t| in denominator
               return (integrator::integrate(fnc_internal_1, xBj, t_abs, 1.0e-5, 1.0e-6).first) / t_abs;
            };
            const double tmp = integrator::integrate(fnc_external, t_min, t_max, 1.0e-5, 1.0e-6).first;
            // Works both if key does and does not exist.
            weight_results[index] += tmp;
            // res += tmp * f._fj(index);
            std::fprintf(stderr, "\r\033[2K"); // Clean line
            std::fprintf(stderr, "Done: %zu, %zu", j, k);
         }
      }

      // Sector 5, both first and second integrals
      for (size_t k = grid.grid_angle._delim_indexes[5]; k < grid.grid_angle._delim_indexes[6]; k++) {
         for (size_t j = 0; j < grid.grid_radius.size; j++) {
            auto r_supp         = grid.grid_radius.get_support_weight_aj(j);
            const double xi_min = std::max(grid.grid_radius._d_info.to_phys_space(r_supp.first), xBj);
            const double xi_max = std::min(grid.grid_radius._d_info.to_phys_space(r_supp.second), 1.0);
            auto index          = f.get_flatten_index(j, k);

            auto fnc_external = [&](double xi) -> double {
               auto fnc_internal_1 = [&](double eta) -> double {
                  auto rhophi = Honeycomb::RnC::from_x123_to_rhophi(xi, -eta - xi, eta);
                  return f._dw_dx3[index](rhophi) / eta;
               };
               auto fnc_internal_2 = [&](double eta) -> double {
                  auto rhophi = Honeycomb::RnC::from_x123_to_rhophi(xi, -eta - xi, eta);
                  return f._dw_dx3[index](rhophi) / (eta + xi);
               };

               // Combines pieces from first and second integral
               return (integrator::integrate(fnc_internal_1, -xi, -xBj, 1.0e-5, 1.0e-6).first +
                       integrator::integrate(fnc_internal_2, -xBj, 0, 1.0e-5, 1.0e-6).first) /
                      xi;
            };
            const double tmp = integrator::integrate(fnc_external, xi_min, xi_max, 1.0e-5, 1.0e-6).first;
            // Works both if key does and does not exist.
            weight_results[index] += tmp;
            // res += tmp * f._fj(index);
            std::fprintf(stderr, "\r\033[2K"); // Clean line
            std::fprintf(stderr, "Done: %zu, %zu", j, k);
         }
      }

      for (const auto &[key, val] : weight_results) {
         res += f._fj(static_cast<long int>(key)) * val;
      }

      const double g2 = res;
      std::fprintf(fp_out, "%.16e\t%+.16e\n", xBj, g2);
      std::fprintf(stderr, "\nDone: %lf, %le\n", xBj, g2);
      x_weight_results.emplace_back(weight_results);
   }
   std::fprintf(stderr, "\n");
   std::fclose(fp_out);
   Honeycomb::timer::mark end = Honeycomb::timer::now();

   Honeycomb::logger(Honeycomb::Logger::INFO,
                     std::format("Full computation took: {:.4e} (ms)", Honeycomb::timer::elapsed_ms(end, begin)));
}