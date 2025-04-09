
#include <honeycomb2/kernel_functions.hpp>
#include <honeycomb2/gauss_kronrod.hpp>

// Local alias for the integrator
using integrator = Honeycomb::Integration::GaussKronrod<Honeycomb::Integration::GK_61>;

namespace Honeycomb
{

double Vplus13::subtracted_integrate(size_t c_a, size_t c_aP, const Grid2D &g)
{
   // Here we are guaranteed to be on the diagonal
   // Retrive external x point
   auto [x1, x2, x3] = g._x123[c_a].v;

   double vmin = 10;
   double vmax = -10;

   double subtr = 0;

   std::vector<long int> w_indexes;

   auto [c_j, c_i] = g.c_get_double_index(c_aP);

   for (size_t c_h = 0; c_h < g.grid_angle._from_ic_to_iw[c_i].size(); c_h++) {
      for (size_t c_l = 0; c_l < g.grid_radius._from_ic_to_iw[c_j].size(); c_l++) {
         long int h = g.grid_angle._from_ic_to_iw[c_i][c_h];
         long int l = g.grid_radius._from_ic_to_iw[c_j][c_l];

         long int aP = g.get_flatten_index(l, h);

         auto [x123_min, x123_max]  = g._x123_minmax[aP];
         auto [x1min, x2min, x3min] = x123_min.v;
         auto [x1max, x2max, x3max] = x123_max.v;

         double loc_vmin = std::max(x1 - x1max, x3min - x3);
         double loc_vmax = std::min(x1 - x1min, x3max - x3);

         vmin = std::min(vmin, loc_vmin);
         vmax = std::max(vmax, loc_vmax);

         subtr += g._w[aP](RnC::from_x123_to_rhophi(x1, x2, x3));

         w_indexes.emplace_back(aP);
      }
   }

   if (x1 == 0.0) {
      return -2.0 * subtr / x3;
   }
   if (x3 == 0.0) {
      return 2.0 * subtr / x1;
   }

   // If support is inverted, then we are outside of it, return 0
   if (vmin >= vmax) return 0;

   if (std::fabs(x2) < 1.0e-14) {

      auto th_x1_mv_integral = [&](double v) -> double {
         double w = 0;
         for (const long int &q : w_indexes)
            w += g._w[q](RnC::from_x123_to_rhophi(x1 - v, 0.0, -x1 + v));

         const double res = w * sq(x1) / sq(sq(v - x1));

         return res;
      };

      double result = 0;

      if (x1 > 0 && vmin <= 0) {
         double upper  = std::min(0.0, vmax);
         result       += integrator::integrate(th_x1_mv_integral, vmin, upper);
      } else if (x1 < 0 && vmax >= 0) {
         double lower  = std::max(vmin, 0.0);
         result       -= integrator::integrate(th_x1_mv_integral, lower, vmax);
      }

      return result;
   }

   const double pref = x1 * x3 / cu(x2);

   // Now we can define the lambdas for the kernel
   auto th_x1_mv_integral = [&](double v) -> double {
      double w = 0;
      for (const long int &q : w_indexes)
         w += g._w[q](RnC::from_x123_to_rhophi(x1 - v, x2, x3 + v));

      w -= subtr;

      const double f = (3 * x1 + x3 - 2 * v) / sq(v - x1);

      return w * f;
   };

   auto th_x3_pv_integral = [&](double v) -> double {
      double w = 0;
      for (const long int &q : w_indexes)
         w += g._w[q](RnC::from_x123_to_rhophi(x1 - v, x2, x3 + v));

      w -= subtr;

      const double f = -(3 * x3 + x1 + 2 * v) / sq(v + x3);

      return w * f;
   };

   double result = 0.0;

   // Now, region by region:
   // First, Theta(x1, -v) = [ x1>=0 && v<=0 ] - [ x1<=0 && v>=0 ]
   // Kernel vanishes for x1=0, so loose inequalities suffice
   if (x1 >= 0 && vmin < 0) {

      const double upper  = std::min(0.0, vmax);
      result             += integrator::integrate(th_x1_mv_integral, vmin, upper);

      if (std::fabs(subtr - 1.0) < 1.0e-12) {
         result +=
             ((upper - vmin) * (x1 + x3)) / ((upper - x1) * (vmin - x1)) - 2 * log(-upper + x1) + 2 * log(-vmin + x1);
      }
   } else if (x1 < 0 && vmax > 0) {

      const double lower  = std::max(0.0, vmin);
      result             -= integrator::integrate(th_x1_mv_integral, lower, vmax); // ! note the - from the Theta

      if (std::fabs(subtr - 1.0) < 1.0e-12) {
         result +=
             -(((lower - vmax) * (x1 + x3)) / ((lower - x1) * (-vmax + x1))) - 2 * log(lower - x1) + 2 * log(vmax - x1);
      }
   }

   if (x3 >= 0 && vmax > 0) {

      const double lower  = std::max(0.0, vmin);
      result             += integrator::integrate(th_x3_pv_integral, lower, vmax);

      if (std::fabs(subtr - 1.0) < 1.0e-12) {
         result +=
             ((lower - vmax) * (x1 + x3)) / ((lower + x3) * (vmax + x3)) + 2 * log(lower + x3) - 2 * log(vmax + x3);
      }

   } else if (x3 < 0 && vmin < 0) {

      const double upper  = std::min(0.0, vmax);
      result             -= integrator::integrate(th_x3_pv_integral, vmin, upper); // ! note the - from the Theta

      if (std::fabs(subtr - 1.0) < 1.0e-12) {
         result += ((upper - vmin) * (x1 + x3)) / ((upper + x3) * (vmin + x3)) + 2 * log((upper + x3) / (vmin + x3));
      }
   }

   return pref * result;
}

double Vplus13::integrate(size_t c_a, size_t aP, const Grid2D &g)
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

   // vmin = max ( x1-x1max, x2min - x2 )
   const double vmin = std::max(x1 - x1max, x3min - x3);
   const double vmax = std::min(x1 - x1min, x3max - x3);

   // If support is inverted, then we are outside of it, return 0
   if (vmin >= vmax) return 0;

   if (std::fabs(x2) < 1.0e-14) {

      auto th_x1_mv_integral = [&](double v) -> double {
         const double w   = g._w[aP](RnC::from_x123_to_rhophi(x1 - v, 0.0, -x1 + v));
         const double res = w * sq(x1) / sq(sq(v - x1));

         return res;
      };

      double result = 0;

      if (x1 > 0 && vmin <= 0) {
         double upper  = std::min(0.0, vmax);
         result       += integrator::integrate(th_x1_mv_integral, vmin, upper);
      } else if (x1 < 0 && vmax >= 0) {
         double lower  = std::max(vmin, 0.0);
         result       -= integrator::integrate(th_x1_mv_integral, lower, vmax);
      }

      return result;

   } else {

      // This is a hack, but result is incorrect for x1=0 or x3=0 otherwise
      if (x1 == 0.0 || x3 == 0.0) {
         return 0.0;
      }

      const double pref = x1 * x3 / cu(x2);

      // Now we can define the lambdas for the kernel
      auto th_x1_mv_integral = [&](double v) -> double {
         const double w = g._w[aP](RnC::from_x123_to_rhophi(x1 - v, x2, x3 + v));
         const double f = (3 * x1 + x3 - 2 * v) / sq(v - x1);

         return w * f;
      };

      auto th_x3_pv_integral = [&](double v) -> double {
         const double w = g._w[aP](RnC::from_x123_to_rhophi(x1 - v, x2, x3 + v));
         const double f = -(3 * x3 + x1 + 2 * v) / sq(v + x3);

         return w * f;
      };

      double result = 0;

      // Now, region by region:
      // First, Theta(x1, -v) = [ x1>=0 && v<=0 ] - [ x1<=0 && v>=0 ]
      // Kernel vanishes for x1=0, so loose inequalities suffice
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

      return pref * result;
   }
}

} // namespace Honeycomb