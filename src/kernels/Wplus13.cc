
#include <honeycomb2/kernel_functions.hpp>
#include <honeycomb2/gauss_kronrod.hpp>

// Local alias for the integrator
using integrator = Honeycomb::Integration::GaussKronrod<Honeycomb::Integration::GK_61>;

namespace Honeycomb
{
double Wplus13::integrate(size_t c_a, size_t aP, const Grid2D &g)
{
   const auto &[x1, x2, x3] = g._x123[c_a].v;

   // Retrieve the support in physical space of weight index aP
   const auto &[x123_min, x123_max]  = g._x123_minmax[aP];
   const auto &[x1min, x2min, x3min] = x123_min.v;
   const auto &[x1max, x2max, x3max] = x123_max.v;

   // x3 is fixed => if outside the support just return 0
   if (x2 < x2min || x2 > x2max) return 0;

   // vmin = max ( x1-x1max, x2min - x2 )
   const double vmin = std::max(x1 - x1max, x3min - x3);
   const double vmax = std::min(x1 - x1min, x3max - x3);

   // If support is inverted, then we are outside of it, return 0
   if (vmin >= vmax) return 0;

   // Now we can define the lambdas for the kernel
   auto th_x1_mv_integral = [&](double v) -> double {
      const double w = g._w[aP](RnC::from_x123_to_rhophi(x1 - v, x2, x3 + v));
      return -0.5 * w;
   };

   auto th_x3_pv_integral = [&](double v) -> double {
      const double w = g._w[aP](RnC::from_x123_to_rhophi(x1 - v, x2, x3 + v));
      return 0.5 * w;
   };

   double result = 0;

   if (x1 > 0 && vmin < 0) {
      const double upper = std::min(0.0, vmax);

      result += integrator::integrate(th_x1_mv_integral, vmin, upper);
   } else if (x1 < 0 && vmax > 0) {
      const double lower = std::max(0.0, vmin);

      result -= integrator::integrate(th_x1_mv_integral, lower, vmax); // ! note the - from the Theta
   }

   if (x3 > 0 && vmax > 0) {
      const double lower = std::max(0.0, vmin);

      result += integrator::integrate(th_x3_pv_integral, lower, vmax);
   } else if (x3 < 0 && vmin < 0) {
      const double upper = std::min(0.0, vmax);

      result -= integrator::integrate(th_x3_pv_integral, vmin, upper); // ! note the - from the Theta
   }

   return result;
}

double Wplus13P23::integrate(size_t c_a, size_t aP, const Grid2D &g)
{
   const auto &[x1, x2, x3] = g._x123[c_a].v;

   // Retrieve the support in physical space of weight index aP
   const auto &[x123_min, x123_max]  = g._x123_minmax[aP];
   const auto &[x1min, x2min, x3min] = x123_min.v;
   const auto &[x1max, x2max, x3max] = x123_max.v;

   // x3 is fixed => if outside the support just return 0
   if (x3 < x2min || x3 > x2max) return 0;

   // vmin = max ( x1-x1max, x2min - x2 )
   const double vmin = std::max(x1 - x1max, x3min - x2);
   const double vmax = std::min(x1 - x1min, x3max - x2);

   // If support is inverted, then we are outside of it, return 0
   if (vmin >= vmax) return 0;

   // Now we can define the lambdas for the kernel
   auto th_x1_mv_integral = [&](double v) -> double {
      const double w = g._w[aP](RnC::from_x123_to_rhophi(x1 - v, x3, x2 + v));
      return -0.5 * w;
   };

   auto th_x2_pv_integral = [&](double v) -> double {
      const double w = g._w[aP](RnC::from_x123_to_rhophi(x1 - v, x3, x2 + v));
      return 0.5 * w;
   };

   double result = 0;

   if (x1 > 0 && vmin < 0) {
      const double upper = std::min(0.0, vmax);

      result += integrator::integrate(th_x1_mv_integral, vmin, upper);
   } else if (x1 < 0 && vmax > 0) {
      const double lower = std::max(0.0, vmin);

      result -= integrator::integrate(th_x1_mv_integral, lower, vmax); // ! note the - from the Theta
   }

   if (x2 > 0 && vmax > 0) {
      const double lower = std::max(0.0, vmin);

      result += integrator::integrate(th_x2_pv_integral, lower, vmax);
   } else if (x2 < 0 && vmin < 0) {
      const double upper = std::min(0.0, vmax);

      result -= integrator::integrate(th_x2_pv_integral, vmin, upper); // ! note the - from the Theta
   }

   return result;
}

} // namespace Honeycomb