
#include <honeycomb2/kernel_functions.hpp>
#include <honeycomb2/gauss_kronrod.hpp>

// Local alias for the integrator
using integrator = Honeycomb::Integration::GaussKronrod<Honeycomb::Integration::GK_61>;

namespace Honeycomb
{
double Hminus12::integrate(size_t c_a, size_t aP, const Grid2D &g)
{
   // Here we are guaranteed to not be on the diagonal

   // Retrive external x point
   const auto &[x1, x2, x3] = g._x123[c_a].v;

   // Retrieve the support in physical space of weight index aP
   const auto &[x123_min, x123_max]  = g._x123_minmax[aP];
   const auto &[x1min, x2min, x3min] = x123_min.v;
   const auto &[x1max, x2max, x3max] = x123_max.v;

   // x3 is fixed => if outside the support just return 0
   if (x3 < x3min || x3 > x3max) return 0;

   const double vmin = std::max(x1 - x2max, x1min - x2);
   const double vmax = std::min(x1 - x2min, x1max - x2);

   // If support is inverted, then we are outside of it, return 0
   if (vmin >= vmax) return 0;

   if (std::fabs(x3) < 1.0e-14) {

      auto th_x1_mv_integral = [&](double v) -> double {
         const double w   = g._w[aP](RnC::from_x123_to_rhophi(x2 + v, -x2 - v, 0.0));
         const double res = w * sq(v) / (2.0 * cu(x1 - v));

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

      // Now we can define the lambdas for the kernel
      auto th_x1_mv_integral = [&](double v) -> double {
         const double w   = g._w[aP](RnC::from_x123_to_rhophi(x2 + v, x1 - v, x3));
         const double res = w * x1 * (2.0 * x2 * (x1 - v) - x1 * (x2 + v)) / (2 * sq(x1 - v) * sq(x3));

         return res;
      };

      auto th_x2_pv_integral = [&](double v) -> double {
         const double w = g._w[aP](RnC::from_x123_to_rhophi(x2 + v, x1 - v, x3));
         return w * sq(x2) / (2 * (x2 + v) * sq(x3));
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

      if (x2 > 0 && vmax > 0) {
         const double lower = std::max(0.0, vmin);

         result += integrator::integrate(th_x2_pv_integral, lower, vmax);
      } else if (x2 < 0 && vmin < 0) {
         const double upper = std::min(0.0, vmax);

         result -= integrator::integrate(th_x2_pv_integral, vmin, upper); // ! note the - from the Theta
      }

      return result;
   }
}

// // This is just to check that the two approaches are equivalent for kernels without subtraction.
// double Hplus12::integrate_v2(size_t c_a, size_t c_aP, const Grid2D &g)
// {
//    // Here we are guaranteed to be on the diagonal
//    // Retrive external x point
//    auto [x1, x2, x3] = g._x123[c_a].v;

//    double vmin = 10;
//    double vmax = -10;

//    std::vector<long int> w_indexes;

//    auto [c_j, c_i] = g.c_get_double_index(c_aP);

//    for (size_t c_h = 0; c_h < g.grid_angle._from_ic_to_iw[c_i].size(); c_h++) {
//       for (size_t c_l = 0; c_l < g.grid_radius._from_ic_to_iw[c_j].size(); c_l++) {
//          long int h = g.grid_angle._from_ic_to_iw[c_i][c_h];
//          long int l = g.grid_radius._from_ic_to_iw[c_j][c_l];

//          long int aP = g.get_flatten_index(l, h);

//          auto [x123_min, x123_max]  = g._x123_minmax[aP];
//          auto [x1min, x2min, x3min] = x123_min.v;
//          auto [x1max, x2max, x3max] = x123_max.v;

//          double loc_vmin = std::max(x1 - x1max, x2min - x2);
//          double loc_vmax = std::min(x1 - x1min, x2max - x2);

//          vmin = std::min(vmin, loc_vmin);
//          vmax = std::max(vmax, loc_vmax);

//          w_indexes.emplace_back(aP);
//       }
//    }

//    // If support is inverted, then we are outside of it, return 0
//    if (vmin >= vmax) return 0;

//    if (std::fabs(x3) < 1.0e-14) {

//       auto th_x1_mv_integral = [&](double v) -> double {
//          double w = 0;
//          for (const long int &q : w_indexes)
//             w += g._w[q](RnC::from_x123_to_rhophi(x1 - v, -x1 + v, 0));

//          const double res = w * v * (v - 2 * x1) / (2.0 * cu(x1 - v));

//          return res;
//       };

//       double result = 0;

//       if (x1 > 0 && vmin <= 0) {
//          double upper  = std::min(0.0, vmax);
//          result       += integrator::integrate(th_x1_mv_integral, vmin, upper);
//       } else if (x1 < 0 && vmax >= 0) {
//          double lower  = std::max(vmin, 0.0);
//          result       -= integrator::integrate(th_x1_mv_integral, lower, vmax);
//       }

//       return result;

//    } else {

//       // Now we can define the lambdas for the kernel
//       auto th_x1_mv_integral = [&](double v) -> double {
//          double w = 0;
//          for (const long int &q : w_indexes)
//             w += g._w[q](RnC::from_x123_to_rhophi(x1 - v, x2 + v, x3));

//          const double res = w * x1 * (x2 - x3) / (2 * (x1 - v) * sq(x3));

//          return res;
//       };

//       auto th_x2_pv_integral = [&](double v) -> double {
//          double w = 0;
//          for (const long int &q : w_indexes)
//             w += g._w[q](RnC::from_x123_to_rhophi(x1 - v, x2 + v, x3));
//          return w * sq(x2) * (x2 - x3 + v) / (2 * sq(x2 + v) * sq(x3));
//       };

//       double result = 0;

//       // Now, region by region:
//       // First, Theta(x1, -v) = [ x1>=0 && v<=0 ] - [ x1<=0 && v>=0 ]
//       // Kernel vanishes for x1=0, so loose inequalities suffice
//       if (x1 > 0 && vmin < 0) {
//          const double upper = std::min(0.0, vmax);

//          result += integrator::integrate(th_x1_mv_integral, vmin, upper);
//       } else if (x1 < 0 && vmax > 0) {
//          const double lower = std::max(0.0, vmin);

//          result -= integrator::integrate(th_x1_mv_integral, lower, vmax); // ! note the - from the Theta
//       }

//       if (x2 > 0 && vmax > 0) {
//          const double lower = std::max(0.0, vmin);

//          result += integrator::integrate(th_x2_pv_integral, lower, vmax);
//       } else if (x2 < 0 && vmin < 0) {
//          const double upper = std::min(0.0, vmax);

//          result -= integrator::integrate(th_x2_pv_integral, vmin, upper); // ! note the - from the Theta
//       }

//       return result;
//    }
// }

} // namespace Honeycomb