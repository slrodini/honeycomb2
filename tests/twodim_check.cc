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

std::vector<double> check_point(const Honeycomb::RnC::Pair &rhophi, const Honeycomb::Discretization &F,
                                const Eigen::VectorXd &fj, model_fnc_t test_fnc, model_fnc_t test_fnc_der,
                                std::FILE *fp)
{
   const Honeycomb::RnC::Triplet x123 = Honeycomb::RnC::from_rhophi_to_x123(rhophi);

   const double exact = test_fnc(x123(0), x123(1), x123(2));

   // const double approx          = F(rhophi);
   // const double approx = F.interpolate_as_weights_v2(rhophi);

   Honeycomb::timer::mark begin  = Honeycomb::timer::now();
   const double approx           = F.interpolate_as_weights_v3(rhophi, fj);
   Honeycomb::timer::mark end    = Honeycomb::timer::now();
   fnc_elapsed                  += Honeycomb::timer::elapsed_ns(end, begin);

   const double dx = 1.0e-7;
   // test_fnc_der(x123(0), x123(1), x123(2))
   const double num_der = (test_fnc(x123(0), -x123(0) - x123(2) - dx, x123(2) + dx) -
                           test_fnc(x123(0), -x123(0) - x123(2) + dx, x123(2) - dx)) /
                          (2.0 * dx);

   const double exact_der = (test_fnc_der != nullptr) ? test_fnc_der(x123(0), x123(1), x123(2)) : num_der;

   begin                   = Honeycomb::timer::now();
   const double inter_der  = F.interpolate_df_dx3_fixed_x1(rhophi, fj);
   end                     = Honeycomb::timer::now();
   fnc_der_elapsed        += Honeycomb::timer::elapsed_ns(end, begin);
   fnc_count++;

   const double err     = std::fabs(approx - exact);
   const double err_der = std::fabs(inter_der - num_der);
   const double err_ex_der =
       std::fabs(exact_der) < 1.0e-10 ? std::fabs(inter_der - exact_der) : std::fabs(1.0 - inter_der / exact_der);

   if (nullptr != fp) {
      fprintf(fp, "%+.16e\t%+.16e\t%+.16e\t%+.16e\t%+.16e\n", rhophi(0), rhophi(1), err, err_der, err_ex_der);
   } else {
      Honeycomb::logger(Honeycomb::Logger::INFO,
                        std::format("{:+.16e}\t{:+.16e}\t{:+.16e}\t{:+.16e}\t{:+.16e}\t{:+.16e}", x123(0), x123(1),
                                    x123(2), err, err_der, err_ex_der));
   }
   return {err, err_der, err_ex_der};
}

int main()
{
   model_fnc_t test     = T_test;
   model_fnc_t test_der = T_test_der;

   const size_t n    = 13;
   const double rmin = 0.01;

   Honeycomb::Grid2D grid = Honeycomb::generate_compliant_Grid2D(n, {rmin, 0.15, 0.65, 1}, {13, 13, 11});

   Honeycomb::Discretization discr(grid);
   Eigen::VectorXd fj = discr(test);

   long int fj_size = discr._grid.c_size_li;

   Honeycomb::logger(Honeycomb::Logger::INFO,
                     std::format("fj size: {:d}, weights size: {:d}", fj_size, discr._grid.size));
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
   const double inter_der = discr.interpolate_df_dx3_fixed_x1(Honeycomb::RnC::from_x123_to_rhophi(x1, x2, x3), fj);
   Honeycomb::logger(Honeycomb::Logger::INFO,
                     std::format("Derivative at ({:+.3f}, {:+.3f}, {:+.3f}) => [E, N, I]: [{:+.8e}, {:+.8e}, {:+.8e}]",
                                 x1, x2, x3, exact_der, num_der, inter_der));

   const double approx = discr.interpolate_as_weights_v3(Honeycomb::RnC::from_x123_to_rhophi(x1, x2, x3), fj);
   const double exact  = test(x1, x2, x3);
   Honeycomb::logger(Honeycomb::Logger::INFO,
                     std::format("Function at   ({:+.3f}, {:+.3f}, {:+.3f}) => [E, I]:    [{:+.10e}, {:+.10e}]", x1, x2,
                                 x3, exact, approx));
   double ave            = 0;
   double ave_der        = 0;
   double ave_ex_der     = 0;
   const size_t nSamples = 1000;
   std::FILE *fp         = fopen("interpolation_checks.dat", "w");

   const double rmax = 0.8;

   for (size_t i = 0; i < nSamples; i++) {
      Honeycomb::RnC::Pair rhophi = {Honeycomb::Random::random_uniform(rmin, rmax),
                                     Honeycomb::Random::random_uniform(0, 6)};

      std::vector<double> cp  = check_point(rhophi, discr, fj, test, test_der, fp);
      const double err        = cp[0];
      const double err_der    = cp[1];
      const double err_ex_der = cp[2];

      const Honeycomb::RnC::Triplet x123 = Honeycomb::RnC::from_rhophi_to_x123(rhophi);

      ave        += err;
      ave_der    += err_der;
      ave_ex_der += err_ex_der;
   }
   Honeycomb::logger(Honeycomb::Logger::INFO,
                     std::format("Average relative error [F, D]: [{:+4e}, {:+.4e}, {:+.4e}]",
                                 ave / static_cast<double>(nSamples), ave_der / static_cast<double>(nSamples),
                                 ave_ex_der / static_cast<double>(nSamples),
                                 ave_ex_der / static_cast<double>(nSamples)));
   fclose(fp);
   Honeycomb::logger(Honeycomb::Logger::INFO, std::format("Average interpolation time:            {:+.6e} ns.",
                                                          fnc_elapsed / static_cast<double>(fnc_count)));

   Honeycomb::logger(Honeycomb::Logger::INFO, std::format("Average derivative interpolation time: {:+.6e} ns.",
                                                          fnc_der_elapsed / static_cast<double>(fnc_count)));
   return 0;
}
