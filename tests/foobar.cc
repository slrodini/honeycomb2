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
double f1_test(double x)
{
   return exp(x);
}

typedef double (*model_fnc_t)(double, double, double);

double fnc_elapsed     = 0;
double fnc_der_elapsed = 0;
size_t fnc_count       = 0;

void check_point(const Honeycomb::RnC::Pair &rhophi, const Honeycomb::Discretization &F, model_fnc_t test_fnc,
                 model_fnc_t test_fnc_der, std::FILE *fp)
{
   const Honeycomb::RnC::Triplet x123 = Honeycomb::RnC::from_rhophi_to_x123(rhophi);

   const double exact = test_fnc(x123(0), x123(1), x123(2));

   // const double approx          = F(rhophi);
   Honeycomb::timer::mark begin = Honeycomb::timer::now();
   const double approx          = F.interpolate_as_weights_v2(rhophi);
   Honeycomb::timer::mark end   = Honeycomb::timer::now();
   fnc_elapsed += Honeycomb::timer::elapsed_ns(end, begin);

   const double dx = 1.0e-7;

   const double num_der = (test_fnc_der != nullptr) ? test_fnc_der(x123(0), x123(1), x123(2))
                                                    : (test_fnc(x123(0), -x123(0) - x123(2) - dx, x123(2) + dx) -
                                                       test_fnc(x123(0), -x123(0) - x123(2) + dx, x123(2) - dx)) /
                                                          (2.0 * dx);

   begin                  = Honeycomb::timer::now();
   const double inter_der = F.interpolate_df_dx3_fixed_x1(rhophi);
   end                    = Honeycomb::timer::now();
   fnc_der_elapsed += Honeycomb::timer::elapsed_ns(end, begin);
   fnc_count++;

   const double err = std::abs(exact) > 1.0e-10 ? std::abs(1 - approx / exact) : std::abs(approx - exact);
   const double err_der =
       std::abs(num_der) > 1.0e-10 ? std::abs(1 - inter_der / num_der) : std::abs(inter_der - num_der);

   if (nullptr != fp) {
      fprintf(fp, "%+.16e\t%+.16e\t%+.16e\t%+.16e\n", rhophi(0), rhophi(1), err, err_der);
   } else {
      Honeycomb::logger(Honeycomb::Logger::INFO, std::format("{:+.16e}\t{:+.16e}\t{:+.16e}\t{:+.16e}\t{:+.16e}",
                                                             x123(0), x123(1), x123(2), err, err_der));
   }
}

void dummy();

int main()
{

   //

   // const double rmax = 0.75;

   std::function<double(double)> r_to_i_s = [](double x) -> double {
      return log(x);
      // const double t = (rmax - rmin) * (x - rmin) / (1 - rmin) + rmin;
      // return log(t);
      // return -log(1 - log(x));
   };
   std::function<double(double)> r_to_i_s_der = [](double x) -> double {
      return 1.0 / x;
      // const double t = (rmax - rmin) * (x - rmin) / (1 - rmin) + rmin;
      // return (rmax - rmin) / (1 - rmin) / t;
      // return 1.0 / (x * (1 - log(x)));
   };

   std::function<double(double)> r_to_p_s = [](double u) -> double {
      return exp(u);
      // const double eu = exp(u);
      // return (eu - rmin) * (1 - rmin) / (rmax - rmin) + rmin;
      // return exp(-expm1(-u));
   };
   std::function<double(double)> r_to_p_s_der = [](double u) -> double {
      return exp(u);
      // const double eu = exp(u);
      // return (eu) * (1 - rmin) / (rmax - rmin);
      // return exp(-expm1(-u) - u);
   };

   const size_t n    = 12;
   const double rmin = 0.01;

   Honeycomb::Grid2D grid = Honeycomb::generate_compliant_Grid2D(n, {rmin, 0.33, 0.66, 1}, {12, 12, 12}, r_to_i_s,
                                                                 r_to_i_s_der, r_to_p_s, r_to_p_s_der);

   model_fnc_t test     = T_test_5;
   model_fnc_t test_der = T_test_5_der;

   Honeycomb::Discretization f(grid, test);
   Honeycomb::logger(Honeycomb::Logger::INFO, std::format("{:d}", f.size_li));
   Honeycomb::logger(
       Honeycomb::Logger::INFO,
       std::format("Estimated kernel sizes in MB: {:f}",
                   static_cast<double>(f.size_li) * static_cast<double>(f.size_li) * 8.0 / (1024.0 * 1024.0)));

   Honeycomb::logger(Honeycomb::Logger::INFO,
                     std::format("Compare with current honeycomb (which needs to store at least the 2 indexes for each "
                                 "kernel as int32_t, plus some "
                                 "additional space for efficient parallelization) MB: {:f}",
                                 503846 * 8.0 * 2 / (1024.0 * 1024.0)));

   const double phi_check                   = grid.grid_angle._coord[15];
   const double rho_check                   = grid.grid_radius._coord[grid.grid_radius._delim_indexes[1]] - 0.001;
   const Honeycomb::RnC::Triplet x123_check = Honeycomb::RnC::from_rhophi_to_x123(rho_check, phi_check);

   // const double x1 = 0.27;
   // const double x3 = 0.52;
   // const double x2 = -x1 - x3;
   const double x1 = x123_check(0);
   const double x3 = x123_check(1);
   const double x2 = x123_check(2);
   Honeycomb::logger(Honeycomb::Logger::INFO, std::format("(x1,x2, x3) = ({:+.16e}, {:+.16e}, {:+.16e})", x1, x2, x3));

   const double dx = 1.0e-5;

   const double num_der   = (test(x1, -x1 - x3 - dx, x3 + dx) - test(x1, -x1 - x3 + dx, x3 - dx)) / (2.0 * dx);
   const double exact_der = test_der == nullptr ? 0 : test_der(x1, -x1 - x3, x3);
   const double inter_der = f.interpolate_df_dx3_fixed_x1(Honeycomb::RnC::from_x123_to_rhophi(x1, x2, x3));
   Honeycomb::logger(
       Honeycomb::Logger::INFO,
       std::format("Derivative at ({:.2f}, {:.2f}, {:.2f}) => NumDer: {:+.10e}\tExactDer: {:+.10e}\tInterDer: {:+.10e}",
                   x1, x2, x3, num_der, exact_der, inter_der));

   const double approx = f.interpolate_as_weights_v2(Honeycomb::RnC::from_x123_to_rhophi(x1, x2, x3));
   const double exact  = test(x1, x2, x3);
   Honeycomb::logger(Honeycomb::Logger::INFO,
                     std::format("Function at ({:.2f}, {:.2f}, {:.2f}) => Exact: {:+.10e}\tInter: {:+.10e}", x1, x2, x3,
                                 exact, approx));

   const size_t nSamples = 1000;
   std::FILE *fp         = fopen("interpolation_checks.dat", "w");
   for (size_t i = 0; i < nSamples; i++) {
      check_point({Honeycomb::Random::random_uniform(rmin, 1), Honeycomb::Random::random_uniform(0, 6)}, f, test,
                  test_der, fp);
   }
   fclose(fp);
   Honeycomb::logger(
       Honeycomb::Logger::INFO,
       std::format("Average interpolation time: {:+.6e} ns; Average derivative interpolation time: {:+.6e} ns.",
                   fnc_elapsed / static_cast<double>(fnc_count), fnc_der_elapsed / static_cast<double>(fnc_count)));
   return 0;
}

void dummy()
{
   const size_t n    = 10;
   const double rmin = 0.01;

   std::function<double(double)> r_to_i_s = [](double x) -> double {
      // return log(x);
      return -acosh(1.0 / cbrt(x));
   };
   std::function<double(double)> r_to_i_s_der = [](double x) -> double {
      // return 1.0 / x;
      const double curt_x  = cbrt(x);
      const double curt_x2 = curt_x * curt_x;
      return 1.0 / (3.0 * x * curt_x * sqrt(1.0 / curt_x2 - 1));
   };

   std::function<double(double)> r_to_p_s = [](double u) -> double {
      // return exp(u);
      const double chu = cosh(-u);
      return 1.0 / (chu * chu * chu);
   };
   std::function<double(double)> r_to_p_s_der = [](double u) -> double {
      // return exp(u);
      const double shu = sinh(-u);
      return 3.0 * shu * shu * shu * tanh(-u);
   };

   Honeycomb::SingleDiscretizationInfo info_radius({rmin, 0.5, 1}, {n, n}, r_to_i_s, r_to_i_s_der, r_to_p_s,
                                                   r_to_p_s_der);
   Honeycomb::SingleDiscretizationInfo info_angle({0, 1, 2, 3, 4, 5, 6}, {n, n, n, n, n, n});

   Honeycomb::Grid2D grid(info_radius, info_angle);
   Honeycomb::Grid grid_radius(info_radius);
   Honeycomb::Grid grid_angle(info_angle);

   Honeycomb::Discretization f(grid, T_test);
   Honeycomb::Discretization1D f_radius(grid_radius, f1_test);
   Honeycomb::Discretization1D f_angle(grid_angle, f1_test);

   // Honeycomb::RnC::Pair rf = Honeycomb::RnC::from_x123_to_rhophi({0.27, -0.52, -0.27 + 0.52});

   // std::FILE *fp = std::fopen("foo.dat", "w");

   // size_t j = 0;
   // for (long int i = 0; i < f.size_li; i++) {
   //    fprintf(fp, "%+.16e\t%+.16e\t%+.16e\t%+.16e\n", f._x123[j](0), f._x123[j](1), f._x123[j](2), f._fj(i));
   //    j++;
   // }
   // fclose(fp);

   Honeycomb::logger(Honeycomb::Logger::INFO,
                     std::format("{:.16e}\t{:.16e}", f({0.27, 0.52, -0.27 - 0.52}), T_test(0.27, 0.52, -0.27 - 0.52)));

   Honeycomb::logger(Honeycomb::Logger::INFO,
                     std::format("Radius only: {:+.16e} at {:.3f}", f_radius(0.59) - f1_test(0.59), 0.59));
   Honeycomb::logger(Honeycomb::Logger::INFO,
                     std::format("Radius only: {:+.16e} at {:.3f}", f_radius(0.29) - f1_test(0.29), 0.29));

   Honeycomb::logger(Honeycomb::Logger::INFO,
                     std::format("Angle only: {:+.16e} at {:.3f}", f_angle(0.29) - f1_test(0.29), 0.29));
   Honeycomb::logger(Honeycomb::Logger::INFO,
                     std::format("Angle only: {:+.16e} at {:.3f}", f_angle(1.29) - f1_test(1.29), 1.29));
   Honeycomb::logger(Honeycomb::Logger::INFO,
                     std::format("Angle only: {:+.16e} at {:.3f}", f_angle(2.29) - f1_test(2.29), 2.29));
   Honeycomb::logger(Honeycomb::Logger::INFO,
                     std::format("Angle only: {:+.16e} at {:.3f}", f_angle(3.29) - f1_test(3.29), 3.29));
   Honeycomb::logger(Honeycomb::Logger::INFO,
                     std::format("Angle only: {:+.16e} at {:.3f}", f_angle(4.29) - f1_test(4.29), 4.29));
   Honeycomb::logger(Honeycomb::Logger::INFO,
                     std::format("Angle only: {:+.16e} at {:.3f}", f_angle(5.29) - f1_test(5.29), 5.29));
}