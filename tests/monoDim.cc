#include "honeycomb2/random_engine.hpp"
#include "honeycomb2/utilities.hpp"
#include <cstdio>
#include <honeycomb2/honeycomb2.hpp>

double test(double x)
{
   return exp(0.5 * (1 - x) * (1 - x));
}

double test_der(double x)
{
   return exp(0.5 * (1 - x) * (1 - x)) * (x - 1);
}

int main()
{
   //
   // Honeycomb::SingleDiscretizationInfo s_info(
   //     {-1.0, 0.0, 1.0}, {16, 16}, false,
   //     [](double x) {
   //        return exp(2 * x);
   //     },
   //     [](double x) {
   //        return 2.0 * exp(2 * x);
   //     },
   //     [](double t) {
   //        return 0.5 * log(t);
   //     },
   //     [](double t) {
   //        return 0.5 / t;
   //     });
   // Honeycomb::Grid grid(s_info);

   // Honeycomb::Grid grid({{-1.0, 0.0, 1.0}, {16, 16}});
   Honeycomb::Grid grid({{-1.0, 1.0}, {16}});
   Honeycomb::Discretization1D discr(grid);

   auto _fj = discr(test);

   std::vector<double> pts = Honeycomb::Random::random_uniform(100, -1.0, 1.0);

   double m_diff     = 0;
   double m_diff_der = 0;
   std::FILE *fp2    = std::fopen("Checks.dat", "w");
   for (auto p : pts) {
      double exact = test(p);
      double inter = discr.interpolate_as_weights(p, _fj);

      double exact_der = test_der(p);
      double inter_der = discr.interpolate_der_as_weights(p, _fj);
      m_diff           = std::max(m_diff, std::fabs(exact - inter));
      m_diff_der       = std::max(m_diff_der, std::fabs(exact_der - inter_der));
      std::fprintf(fp2, "%.16e\t%.16e\t%.16e\t%.16e\t%.16e\n", p, exact - inter, exact_der - inter_der,
                   exact_der, inter_der);
   }
   std::fclose(fp2);
   Honeycomb::logger(Honeycomb::Logger::INFO, std::format("{:.16e}", m_diff));
   Honeycomb::logger(Honeycomb::Logger::INFO, std::format("{:.16e}", m_diff_der));

   return 0;
}