#include <honeycomb2/alpha_s.hpp>
#include <honeycomb2/utilities.hpp>

namespace Honeycomb
{

std::function<double(double)> GetAlphaS_o_4pi(std::array<double, 6> thresholds, double Q2ref,
                                              double alpha_s_ref, double Nc)
{
   size_t nf_min = 0, nf_ref = 0;

   while (nf_min < 6 && is_near(thresholds[nf_min], 0.0)) {
      nf_min++;
   }
   while (nf_ref < 6 && Q2ref > thresholds[nf_ref]) {
      nf_ref++;
   }
   logger(Logger::INFO,
          std::format("Alpha_s: min value of nf: {:d}; nf at ref scale: {:d}", nf_min, nf_ref));
   std::map<int, double> beta_0;
   for (size_t nf = 0; nf <= 6; nf++) {
      beta_0[nf] = 11.0 * Nc / 3.0 - 2.0 * static_cast<double>(nf) / 3.0;
   }

   std::array<double, 6> t_val;
   for (size_t i = 0; i < 6; i++) {
      if (is_near(thresholds[i], 0.0)) {
         t_val[i] = -std::numeric_limits<double>::infinity();
      } else {
         t_val[i] = log(thresholds[i]);
      }
   }

   // as = (beta_0 (t - L))^(-1), L = Log(Lambda^2)
   std::map<int, double> L;
   L[static_cast<int>(nf_ref)] = log(Q2ref) - 4.0 * M_PI / (alpha_s_ref * beta_0[nf_ref]);
   for (int nf = static_cast<int>(nf_ref) - 1; nf >= static_cast<int>(nf_min); nf--) {
      double tmp = t_val[nf] + (L[nf + 1] - t_val[nf]) * beta_0[nf + 1] / beta_0[nf];
      L[nf]      = tmp;
   }
   for (int nf = nf_ref; nf < 6; nf++) {
      L[nf + 1] = t_val[nf] + (L[nf] - t_val[nf]) * beta_0[nf] / beta_0[nf + 1];
   }

   auto as = [L, t_val, beta_0](double t) -> double {
      int nf = 0;
      while (nf < 6 && t >= t_val[nf]) {
         nf++;
      }
      return 1.0 / (beta_0.at(nf) * (t - L.at(nf)));
   };
   return as;
}
} // namespace Honeycomb