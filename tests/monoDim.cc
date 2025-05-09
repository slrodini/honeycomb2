#include "honeycomb2/discretization.hpp"
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

double test_2(double x)
{
   return exp(2.0 * pow(0.01, (1 - x) * 0.5));
}

double test_2_der(double x)
{
   const double a = 0.01;
   return -exp(2.0 * pow(a, (1 - x) * 0.5)) * pow(a, (1 - x) * 0.5) * log(a);
}

double test_3(double x)
{
   double a = 19;
   return pow(1 - x, a);
}

double test_3_der(double x)
{
   double a = 19;
   return pow(1 - x, a - 1) * (-a);
}

typedef double (*model_t)(double);

void populate_grids(std::vector<Honeycomb::Grid> &grids, std::vector<Honeycomb::Discretization1D> &discrs,
                    double xmin, double xmax);

int main(int argc, char **argv)
{

   model_t f;
   model_t df;
   if (argc <= 1) {
      f  = test;
      df = test_der;
   } else {
      int i = std::atoi(argv[1]);
      switch (i) {
      case 2:
         f  = test_2;
         df = test_2_der;
         break;
      case 3:
         f  = test_3;
         df = test_3_der;
         break;
      case 1:
      default:
         f  = test;
         df = test_der;
         break;
      }
   }
   //

   double xmin = 0.001;
   double xmax = 0.999;

   size_t n_pts = 15;

   Honeycomb::SingleDiscretizationInfo s_info_1(
       {xmin, (xmin + xmax) * 0.5, xmax}, {n_pts, n_pts}, false,
       //  {xmin, xmax}, {n_pts}, false,
       [](double x) {
          return -log(1 - x);
       },
       [](double x) {
          return 1.0 / (1.0 - x);
       },
       [](double t) {
          return 1.0 - exp(-t);
       },
       [](double t) {
          return exp(-t);
       });

   Honeycomb::SingleDiscretizationInfo s_info_2(
       {xmin, (xmin + xmax) * 0.5, xmax}, {n_pts, n_pts}, false,
       //  {xmin, xmax}, {n_pts}, false,
       [](double x) {
          return log(x);
       },
       [](double x) {
          return 1.0 / x;
       },
       [](double t) {
          return exp(t);
       },
       [](double t) {
          return exp(t);
       });

   Honeycomb::Grid grids_0(s_info_1);
   Honeycomb::Grid grids_1(s_info_2);
   // Honeycomb::Grid grids_2({{xmin, 0.5 * (xmin + xmax), xmax}, {n_pts, n_pts}});
   Honeycomb::Grid grids_2({{xmin, xmax}, {n_pts}});

   Honeycomb::Discretization1D discrs_0(grids_0);
   Honeycomb::Discretization1D discrs_1(grids_1);
   Honeycomb::Discretization1D discrs_2(grids_2);

   Eigen::VectorXd _fj_0 = discrs_0(f);
   Eigen::VectorXd _fj_1 = discrs_1(f);
   Eigen::VectorXd _fj_2 = discrs_2(f);

   std::vector<double> pts = Honeycomb::Random::random_uniform(100, xmin, xmax);
   pts.push_back(xmax);
   // std::vector<double> pts = grids_0._coord;

   double m_diff_0     = 0;
   double m_diff_der_0 = 0;

   double m_diff_1     = 0;
   double m_diff_der_1 = 0;

   double m_diff_2     = 0;
   double m_diff_der_2 = 0;

   double m_diff_rel_0     = 0;
   double m_diff_rel_der_0 = 0;

   double m_diff_rel_1     = 0;
   double m_diff_rel_der_1 = 0;

   double m_diff_rel_2     = 0;
   double m_diff_rel_der_2 = 0;
   std::FILE *fp2          = std::fopen("Checks.dat", "w");
   for (auto p : pts) {
      std::fprintf(fp2, "%.16e\t", p);

      double exact     = f(p);
      double exact_der = df(p);
      double inter, inter_der;

      inter        = discrs_0.interpolate_as_weights(p, _fj_0);
      inter_der    = discrs_0.interpolate_der_as_weights(p, _fj_0);
      m_diff_0     = std::max(m_diff_0, std::fabs(exact - inter));
      m_diff_der_0 = std::max(m_diff_der_0, std::fabs(exact_der - inter_der));
      m_diff_rel_0 = std::max(m_diff_rel_0, fabs(exact) < 1.0e-12 ? std::fabs(exact - inter)
                                                                  : std::fabs(1.0 - inter / exact));
      m_diff_rel_der_0
          = std::max(m_diff_rel_der_0, fabs(exact_der) < 1.0e-12 ? std::fabs(exact_der - inter_der)
                                                                 : std::fabs(1.0 - inter_der / exact_der));
      std::fprintf(fp2, "%.16e\t%.16e\t", exact_der, inter_der);

      inter        = discrs_1.interpolate_as_weights(p, _fj_1);
      inter_der    = discrs_1.interpolate_der_as_weights(p, _fj_1);
      m_diff_1     = std::max(m_diff_1, std::fabs(exact - inter));
      m_diff_der_1 = std::max(m_diff_der_1, std::fabs(exact_der - inter_der));
      m_diff_rel_1 = std::max(m_diff_rel_1, fabs(exact) < 1.0e-12 ? std::fabs(exact - inter)
                                                                  : std::fabs(1.0 - inter / exact));
      m_diff_rel_der_1
          = std::max(m_diff_rel_der_1, fabs(exact_der) < 1.0e-12 ? std::fabs(exact_der - inter_der)
                                                                 : std::fabs(1.0 - inter_der / exact_der));
      std::fprintf(fp2, "%.16e\t%.16e\t", exact_der, inter_der);

      inter        = discrs_2.interpolate_as_weights(p, _fj_2);
      inter_der    = discrs_2.interpolate_der_as_weights(p, _fj_2);
      m_diff_2     = std::max(m_diff_2, std::fabs(exact - inter));
      m_diff_der_2 = std::max(m_diff_der_2, std::fabs(exact_der - inter_der));
      m_diff_rel_2 = std::max(m_diff_rel_2, fabs(exact) < 1.0e-12 ? std::fabs(exact - inter)
                                                                  : std::fabs(1.0 - inter / exact));
      m_diff_rel_der_2
          = std::max(m_diff_rel_der_2, fabs(exact_der) < 1.0e-12 ? std::fabs(exact_der - inter_der)
                                                                 : std::fabs(1.0 - inter_der / exact_der));
      std::fprintf(fp2, "%.16e\t%.16e\t", exact_der, inter_der);

      std::fprintf(fp2, "\n");
   }
   std::fclose(fp2);
   Honeycomb::logger(Honeycomb::Logger::INFO, std::format("{:.16e}", m_diff_0));
   Honeycomb::logger(Honeycomb::Logger::INFO, std::format("{:.16e}", m_diff_der_0));
   Honeycomb::logger(Honeycomb::Logger::INFO, std::format("{:.16e}", m_diff_rel_0));
   Honeycomb::logger(Honeycomb::Logger::INFO, std::format("{:.16e}", m_diff_rel_der_0));
   Honeycomb::logger(Honeycomb::Logger::INFO, "=================================================");

   Honeycomb::logger(Honeycomb::Logger::INFO, std::format("{:.16e}", m_diff_1));
   Honeycomb::logger(Honeycomb::Logger::INFO, std::format("{:.16e}", m_diff_der_1));
   Honeycomb::logger(Honeycomb::Logger::INFO, std::format("{:.16e}", m_diff_rel_1));
   Honeycomb::logger(Honeycomb::Logger::INFO, std::format("{:.16e}", m_diff_rel_der_1));
   Honeycomb::logger(Honeycomb::Logger::INFO, "=================================================");

   Honeycomb::logger(Honeycomb::Logger::INFO, std::format("{:.16e}", m_diff_2));
   Honeycomb::logger(Honeycomb::Logger::INFO, std::format("{:.16e}", m_diff_der_2));
   Honeycomb::logger(Honeycomb::Logger::INFO, std::format("{:.16e}", m_diff_rel_2));
   Honeycomb::logger(Honeycomb::Logger::INFO, std::format("{:.16e}", m_diff_rel_der_2));
   Honeycomb::logger(Honeycomb::Logger::INFO, "=================================================");

   return 0;
}

void populate_grids(std::vector<Honeycomb::Grid> &grids, std::vector<Honeycomb::Discretization1D> &discrs,
                    double xmin, double xmax)
{
   Honeycomb::SingleDiscretizationInfo s_info_1(
       {xmin, xmax}, {15}, false,
       [](double x) {
          return exp(2 * x);
       },
       [](double x) {
          return 2.0 * exp(2 * x);
       },
       [](double t) {
          return 0.5 * log(t);
       },
       [](double t) {
          return 0.5 / t;
       });

   Honeycomb::SingleDiscretizationInfo s_info_2(
       {xmin, xmax}, {15}, false,
       [](double x) {
          return log(x);
       },
       [](double x) {
          return 1.0 / x;
       },
       [](double t) {
          return exp(t);
       },
       [](double t) {
          return exp(t);
       });

   grids.push_back(Honeycomb::Grid(s_info_1));
   discrs.push_back(Honeycomb::Discretization1D(grids[0]));

   grids.push_back(Honeycomb::Grid(s_info_2));
   discrs.push_back(Honeycomb::Discretization1D(grids[1]));

   grids.push_back(Honeycomb::Grid({{xmin, xmax}, {16}}));
   discrs.push_back(Honeycomb::Discretization1D(grids[2]));
   return;
}
