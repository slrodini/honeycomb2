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

double T_test_2(double x1, double x2, double x3)
{
   return (1.0 - cos(T_test(x1, x2, x3))) * 2 * sin(M_PI * x2) /
          std::max(std::max(std::fabs(x1), std::fabs(x2)), std::fabs(x3));
}

double T_test_3(double x1, double x2, double x3)
{
   return (1.0 - cos(T_test(x1, x2, x3))) * 2 * sin(M_PI * x2) / (x1 * x1 + x2 * x2 + x3 * x3);
}

double T_test_3_der(double x1, double x2, double x3)
{
   return (-2 * M_PI * cos(M_PI * x2) * (1 - cos((1 - pow(x1, 2)) * (1 - pow(x2, 2)) * (1 - pow(x3, 2))))) /
              (pow(x1, 2) + pow(x2, 2) + pow(x3, 2)) -
          (2 * (-2 * x2 + 2 * x3) * (1 - cos((1 - pow(x1, 2)) * (1 - pow(x2, 2)) * (1 - pow(x3, 2)))) *
           sin(M_PI * x2)) /
              pow(pow(x1, 2) + pow(x2, 2) + pow(x3, 2), 2) +
          (2 * (-2 * (1 - pow(x1, 2)) * (1 - pow(x2, 2)) * x3 + 2 * (1 - pow(x1, 2)) * x2 * (1 - pow(x3, 2))) *
           sin(M_PI * x2) * sin((1 - pow(x1, 2)) * (1 - pow(x2, 2)) * (1 - pow(x3, 2)))) /
              (pow(x1, 2) + pow(x2, 2) + pow(x3, 2));
}

double T_test_4(double x1, double x2, double x3)
{
   return (1.0 - cos(T_test(x1, x2, x3))) * 2 * sin(M_PI * x2);
}

double T_test_4_der(double x1, double x2, double x3)
{
   return -2 * Pi * Cos(Pi * x2) * (1 - Cos((1 - Power(x1, 2)) * (1 - Power(x2, 2)) * (1 - Power(x3, 2)))) +
          2 * (-2 * (1 - Power(x1, 2)) * (1 - Power(x2, 2)) * x3 + 2 * (1 - Power(x1, 2)) * x2 * (1 - Power(x3, 2))) *
              Sin(Pi * x2) * Sin((1 - Power(x1, 2)) * (1 - Power(x2, 2)) * (1 - Power(x3, 2)));
}

double T_test_5(double x1, double x2, double x3)
{
   return T_test(x1, x2, x3) * (log(x1 * x1 + x2 * x2 + x3 * x3) + cos(M_PI * x1 * x3));
}

double T_test_5_der(double x1, double x2, double x3)
{
   return -2 * (1 - Power(x1, 2)) * (1 - Power(x2, 2)) * x3 *
              (Cos(Pi * x1 * x3) + Log(Power(x1, 2) + Power(x2, 2) + Power(x3, 2))) +
          2 * (1 - Power(x1, 2)) * x2 * (1 - Power(x3, 2)) *
              (Cos(Pi * x1 * x3) + Log(Power(x1, 2) + Power(x2, 2) + Power(x3, 2))) +
          (1 - Power(x1, 2)) * (1 - Power(x2, 2)) * (1 - Power(x3, 2)) *
              ((-2 * x2 + 2 * x3) / (Power(x1, 2) + Power(x2, 2) + Power(x3, 2)) - Pi * x1 * Sin(Pi * x1 * x3));
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

double f1_test(double x)
{
   return exp(x);
}

typedef double (*model_fnc_t)(double, double, double);

double fnc_elapsed     = 0;
double fnc_der_elapsed = 0;
size_t fnc_count       = 0;

std::pair<double, double> check_point(const Honeycomb::RnC::Pair &rhophi, const Honeycomb::Discretization &F,
                                      model_fnc_t test_fnc, model_fnc_t test_fnc_der, std::FILE *fp)
{
   const Honeycomb::RnC::Triplet x123 = Honeycomb::RnC::from_rhophi_to_x123(rhophi);

   const double exact = test_fnc(x123(0), x123(1), x123(2));

   // const double approx          = F(rhophi);
   // const double approx = F.interpolate_as_weights_v2(rhophi);

   Honeycomb::timer::mark begin  = Honeycomb::timer::now();
   const double approx           = F.interpolate_as_weights_v3(rhophi);
   Honeycomb::timer::mark end    = Honeycomb::timer::now();
   fnc_elapsed                  += Honeycomb::timer::elapsed_ns(end, begin);

   const double dx = 1.0e-7;

   const double num_der = (test_fnc_der != nullptr) ? test_fnc_der(x123(0), x123(1), x123(2))
                                                    : (test_fnc(x123(0), -x123(0) - x123(2) - dx, x123(2) + dx) -
                                                       test_fnc(x123(0), -x123(0) - x123(2) + dx, x123(2) - dx)) /
                                                          (2.0 * dx);

   begin                   = Honeycomb::timer::now();
   const double inter_der  = F.interpolate_df_dx3_fixed_x1(rhophi);
   end                     = Honeycomb::timer::now();
   fnc_der_elapsed        += Honeycomb::timer::elapsed_ns(end, begin);
   fnc_count++;

   // const double err = std::abs(exact) > 1.0e-10 ? std::abs(1 - approx / exact) : std::abs(approx - exact);
   // const double err_der =
   //     std::abs(num_der) > 1.0e-10 ? std::abs(1 - inter_der / num_der) : std::abs(inter_der - num_der);

   const double err     = std::abs(approx - exact);
   const double err_der = std::abs(inter_der - num_der);

   if (nullptr != fp) {
      fprintf(fp, "%+.16e\t%+.16e\t%+.16e\t%+.16e\n", rhophi(0), rhophi(1), err, err_der);
   } else {
      Honeycomb::logger(Honeycomb::Logger::INFO, std::format("{:+.16e}\t{:+.16e}\t{:+.16e}\t{:+.16e}\t{:+.16e}",
                                                             x123(0), x123(1), x123(2), err, err_der));
   }
   return {err, err_der};
}

int main()
{
   model_fnc_t test     = T_test_6;
   model_fnc_t test_der = T_test_6_der;

   const size_t n    = 13;
   const double rmin = 0.01;

   Honeycomb::Grid2D grid = Honeycomb::generate_compliant_Grid2D(n, {rmin, 0.15, 0.65, 1}, {13, 13, 11});

   Honeycomb::Discretization f(grid, test);
   long int fj_size = f._grid.c_size_li;

   Honeycomb::logger(Honeycomb::Logger::INFO, std::format("fj size: {:d}, weights size: {:d}", fj_size, f._grid.size));
   Honeycomb::logger(Honeycomb::Logger::INFO, std::format("Estimated kernel sizes in MB:      {:f}",
                                                          static_cast<double>(fj_size) * static_cast<double>(fj_size) *
                                                              8.0 / (1024.0 * 1024.0)));

   //
   // Current honeycomb needs to store at least the 2 indexes for each kernel as int32_t, plus some additional space for
   // efficient parallelization
   //
   Honeycomb::logger(Honeycomb::Logger::INFO,
                     std::format("Compare with current honeycomb MB: {:f}", 503846 * 8.0 * 2 / (1024.0 * 1024.0)));

   const double phi_check                   = grid.grid_angle._coord[15];
   const double rho_check                   = grid.grid_radius._coord[grid.grid_radius._delim_indexes[1]] - 0.001;
   const Honeycomb::RnC::Triplet x123_check = Honeycomb::RnC::from_rhophi_to_x123(rho_check, phi_check);

   const double x1 = x123_check(0);
   const double x3 = x123_check(1);
   const double x2 = x123_check(2);

   const double dx = 1.0e-5;

   const double num_der   = (test(x1, -x1 - x3 - dx, x3 + dx) - test(x1, -x1 - x3 + dx, x3 - dx)) / (2.0 * dx);
   const double exact_der = test_der == nullptr ? 0 : test_der(x1, -x1 - x3, x3);
   const double inter_der = f.interpolate_df_dx3_fixed_x1(Honeycomb::RnC::from_x123_to_rhophi(x1, x2, x3));
   Honeycomb::logger(Honeycomb::Logger::INFO,
                     std::format("Derivative at ({:+.3f}, {:+.3f}, {:+.3f}) => [E, N, I]: [{:+.8e}, {:+.8e}, {:+.8e}]",
                                 x1, x2, x3, exact_der, num_der, inter_der));

   const double approx = f.interpolate_as_weights_v3(Honeycomb::RnC::from_x123_to_rhophi(x1, x2, x3));
   const double exact  = test(x1, x2, x3);
   Honeycomb::logger(Honeycomb::Logger::INFO,
                     std::format("Function at   ({:+.3f}, {:+.3f}, {:+.3f}) => [E, I]:    [{:+.10e}, {:+.10e}]", x1, x2,
                                 x3, exact, approx));
   double ave            = 0;
   double ave_der        = 0;
   const size_t nSamples = 1000;
   std::FILE *fp         = fopen("interpolation_checks.dat", "w");
   for (size_t i = 0; i < nSamples; i++) {
      auto [err, err_der] = check_point(
          {Honeycomb::Random::random_uniform(rmin, 1), Honeycomb::Random::random_uniform(0, 6)}, f, test, test_der, fp);
      ave     += err;
      ave_der += err_der;
   }
   Honeycomb::logger(Honeycomb::Logger::INFO,
                     std::format("Average relative error [F, D]: [{:+4e}, {:+.4e}]",
                                 ave / static_cast<double>(nSamples), ave_der / static_cast<double>(nSamples)));
   fclose(fp);
   Honeycomb::logger(Honeycomb::Logger::INFO, std::format("Average interpolation time:            {:+.6e} ns.",
                                                          fnc_elapsed / static_cast<double>(fnc_count)));

   Honeycomb::logger(Honeycomb::Logger::INFO, std::format("Average derivative interpolation time: {:+.6e} ns.",
                                                          fnc_der_elapsed / static_cast<double>(fnc_count)));
   return 0;
}
