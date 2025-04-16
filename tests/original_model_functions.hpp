
#include <cmath>

namespace orignal_models
{
inline double max3(double a, double b, double c)
{
   return std::max(a, std::max(b, c));
}

double model_zero_function(double, double, double)
{
   return 0.0;
}

double Tu_test(double x1, double x2, double x3)
{

   return (1 - x1 * x1) * (1 - x2 * x2) * (1 - x3 * x3) * cos(4.0 * x2);
}
double Td_test(double x1, double x2, double x3)
{

   double temp = (1 - x1 * x1) * (1 - x2 * x2) * (1 - x3 * x3);
   return (2 - cos(3 * M_PI * temp)) * temp;
}
double Ts_test(double x1, double x2, double x3)
{

   return -0.3 * Td_test(x1, x2, x3);
}

double DTu_test(double x1, double x2, double x3)
{

   return (sin(x2 * M_PI) + 4 * (x1 * x1 - x3 * x3)) * (1 - x1 * x1) * (1 - x2 * x2) * (1 - x3 * x3);
}
double DTd_test(double x1, double x2, double x3)
{

   double r = max3(fabs(x1), fabs(x2), fabs(x3));
   return sin(x2 * M_PI) * (2 - 2 * cos(3 * M_PI * (1 - x1 * x1) * (1 - x2 * x2) * (1 - x3 * x3))) / sqrt(r);
}
double DTs_test(double x1, double x2, double x3)
{

   return -0.3 * DTd_test(x1, x2, x3);
}

double Hu_test(double x1, double x2, double x3)
{

   double r = max3(fabs(x1), fabs(x2), fabs(x3));
   return sin(x2 * M_PI) * (1 - x1 * x1) * (1 - x2 * x2) * (1 - x3 * x3) / (sqrt(r));
}
double Hd_test(double x1, double x2, double x3)
{
   (void)x1;
   (void)x2;
   (void)x3;

   return 0;
}
double Hs_test(double x1, double x2, double x3)
{
   (void)x1;
   (void)x2;
   (void)x3;

   return 0;
}

double Eu_test(double x1, double x2, double x3)
{

   return (1 - x1 * x1) * (1 - x2 * x2) * (1 - x3 * x3);
}
double Ed_test(double x1, double x2, double x3)
{
   (void)x1;
   (void)x2;
   (void)x3;

   return 0;
}
double Es_test(double x1, double x2, double x3)
{
   (void)x1;
   (void)x2;
   (void)x3;

   return 0;
}

double TFp_test(double x1, double x2, double x3)
{

   double r = max3(fabs(x1), fabs(x2), fabs(x3));
   return sin(x1 - x3) * (1 - x1 * x1) * (1 - x2 * x2) * (1 - x3 * x3) * sqrt(r);
}
double TFm_test(double x1, double x2, double x3)
{

   double r = max3(fabs(x1), fabs(x2), fabs(x3));
   return cos(x1 - x3) * (1 - x1 * x1) * (1 - x2 * x2) * (1 - x3 * x3) * sqrt(r);
}
} // namespace orignal_models