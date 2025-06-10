#include <honeycomb2/honeycomb2.hpp>
#include <honeycomb2/honeycomb2_c_api.h>
#include <test_config.hpp>

int check_initial_condition();

int main()
{
   std::string file = TESTS_PATH "/example.config";

   hc2_fi_set_up_(file.c_str(), file.size());

   hc2_fi_evolve_();
   hc2_fi_unload_();

   return 0;
}
