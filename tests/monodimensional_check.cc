#include <honeycomb2/honeycomb2.hpp>

double f1_test(double x)
{
   return exp(x);
}

double f1_test_per(double x)
{
   return sin(x);
}

typedef double (*model_fnc_t)(double);

double fnc_elapsed     = 0;
double fnc_der_elapsed = 0;
size_t fnc_count       = 0;

double check_point(const double &x, const Honeycomb::Discretization1D &F, model_fnc_t test_fnc, std::FILE *fp)
{
   const double exact = test_fnc(x);

   // const double approx          = F(rhophi);
   Honeycomb::timer::mark begin  = Honeycomb::timer::now();
   const double approx           = F.interpolate_as_weights(x);
   Honeycomb::timer::mark end    = Honeycomb::timer::now();
   fnc_elapsed                  += Honeycomb::timer::elapsed_ns(end, begin);

   fnc_count++;

   const double err = std::abs(exact) > 1.0e-10 ? std::abs(1 - approx / exact) : std::abs(approx - exact);

   if (nullptr != fp) {
      fprintf(fp, "%+.16e\t%+.16e\n", x, err);
   } else {
      Honeycomb::logger(Honeycomb::Logger::INFO, std::format("{:+.16e}\t{:+.16e}", x, err));
   }
   return err;
}

int main()
{

   const size_t n    = 16;
   const double rmin = 0.01;

   Honeycomb::Grid grid({{rmin, 0.5, 1}, {n, n}});
   Honeycomb::Grid grid_per({{0, M_PI, 2 * M_PI}, {n, n}, true});

   model_fnc_t test     = f1_test;
   model_fnc_t test_per = f1_test_per;

   Honeycomb::Discretization1D f(grid, test);
   Honeycomb::Discretization1D f_per(grid_per, test_per);

   Honeycomb::logger(Honeycomb::Logger::INFO, std::format("Grid sizes: {:d}, {:d}; fj sizes: {:d}, {:d}", f.size,
                                                          f_per.size, f._fj.size(), f_per._fj.size()));

   // for (size_t i = 0; i < f_per._grid.size; i++) {
   //    Honeycomb::logger(Honeycomb::Logger::INFO, std::format("{:.16e}", f_per._fj(grid_per._from_iw_to_ic[i])));
   // }
   // exit(0);

   double ave            = 0;
   double ave_per        = 0;
   const size_t nSamples = 1000;
   std::FILE *fp         = fopen("interpolation_checks.dat", "w");
   std::FILE *fp_per     = fopen("interpolation_checks_per.dat", "w");
   for (size_t i = 0; i < nSamples; i++) {
      ave     += check_point(Honeycomb::Random::random_uniform(rmin, 1), f, test, fp);
      ave_per += check_point(Honeycomb::Random::random_uniform(0, 2 * M_PI), f_per, test_per, fp_per);
   }
   std::fclose(fp);
   std::fclose(fp_per);
   Honeycomb::logger(Honeycomb::Logger::INFO,
                     std::format("Average relative error. Exp(x): {:+16e}\t Sin(x) with periodic grid: {:+.16e}",
                                 ave / static_cast<double>(nSamples), ave_per / static_cast<double>(nSamples)));
   Honeycomb::logger(
       Honeycomb::Logger::INFO,
       std::format("Average interpolation time: {:+.6e} ns; Average derivative interpolation time: {:+.6e} ns.",
                   fnc_elapsed / static_cast<double>(fnc_count), fnc_der_elapsed / static_cast<double>(fnc_count)));
   return 0;
}
