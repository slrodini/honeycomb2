#include <honeycomb2/honeycomb2.hpp>

int main()
{
   std::array<double, 6> thresholds = {0, 0, 0, 1.6129, 17.4724, 1.0e+6};

   auto as       = Honeycomb::GetAlphaS_o_4pi(thresholds);
   std::FILE *fp = std::fopen("alpha_s.dat", "w");
   for (double t = 0.01; t < 2 * log(thresholds[5]); t *= 1.05) {
      std::fprintf(fp, "%.16e\t%.16e\n", t, as(t));
   }
   std::fclose(fp);
   return 0;
}