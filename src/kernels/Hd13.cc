
#include <honeycomb2/kernel_functions.hpp>
#include <honeycomb2/gauss_kronrod.hpp>

// Local alias for the integrator
using integrator = Honeycomb::Integration::GaussKronrod<Honeycomb::Integration::GK_61>;

namespace Honeycomb
{
double Hd13::integrate(size_t c_a, size_t aP, const Grid2D &g)
{
   // Here we are guaranteed to not be on the diagonal

   // Retrive external x point
   const auto &[x1, x2, x3] = g._x123[c_a].v;

   // Retrieve the support in physical space of weight index aP
   const auto &[x123_min, x123_max]  = g._x123_minmax[aP];
   const auto &[x1min, x2min, x3min] = x123_min.v;
   const auto &[x1max, x2max, x3max] = x123_max.v;

   // x3 is fixed => if outside the support just return 0
   if (x2 < x2min || x2 > x2max) return 0;
   if (x2 == 0.0) return 0.0;

   // vmin = max ( x1-x1max, x2min - x2 )
   const double vmin = std::max(x1 - x1max, x3min - x3);
   const double vmax = std::min(x1 - x1min, x3max - x3);

   // If support is inverted, then we are outside of it, return 0
   if (vmin >= vmax) return 0;

   const double m_theta = (x1 > 0 && x3 > 0) ? -1.0 : ((x1 < 0 && x3 < 0) ? +1.0 : 0.0);
   if (m_theta == 0.0) return 0.0;

   auto th_x1_mv_integral = [&](double v) -> double {
      const double w = g._w[aP](RnC::from_x123_to_rhophi(x1 - v, x2, x3 + v));
      return w;
   };

   const double res = integrator::integrate(th_x1_mv_integral, vmin, vmax);
   return m_theta * x1 * x3 * res / sq(x2);
}

} // namespace Honeycomb