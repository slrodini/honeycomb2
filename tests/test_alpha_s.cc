#include "utilities.hpp"
#include <honeycomb2/honeycomb2.hpp>
#include "test_rd_alpha_s.hpp"

#define GEN_REGRESSION 0

void generate_regression(std::function<double(double)> &as,
                         const std::array<double, 6> &thresholds);

int main()
{
   std::array<double, 6> thresholds = {0, 0, 0, 1.6129, 17.4724, 1.0e+6};

   auto as = Honeycomb::GetAlphaS_o_4pi(thresholds);
   if (GEN_REGRESSION == 1) {
      generate_regression(as, thresholds);
      return 0;
   }

   int ret_val = 0;
   for (const auto [t, as_r] : as_reference) {
      if (!Honeycomb::is_near(as_r, as(t))) {
         ret_val += 1;
         Honeycomb::logger(Honeycomb::Logger::WARNING, std::format("Offending point: t={:.12e}, "
                                                                   "as_ref={:.12e}, "
                                                                   "as={:.12e}",
                                                                   t, as_r, as(t)));
      }
   }
   if (ret_val > 0 && ret_val % 255 == 0) {
      ret_val = -1;
   }
   return ret_val;
}

void generate_regression(std::function<double(double)> &as, const std::array<double, 6> &thresholds)
{
   std::FILE *fp = std::fopen("alpha_s.hpp", "w");
   std::fprintf(fp, "#include <vector>\n#include<utility>\n\n");
   std::fprintf(fp, "std::vector<std::pair<double, double>> as_reference = {\n");
   for (double t = 0.01; t < 2 * log(thresholds[5]); t *= 1.05) {
      std::fprintf(fp, "{%.16e, %.16e},\n", t, as(t));
   }
   std::fprintf(fp, "};\n");

   std::fclose(fp);
}
