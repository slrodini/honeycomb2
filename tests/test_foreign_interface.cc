#include <Honeycomb2/honeycomb2.hpp>
#include <Honeycomb2/honeycomb2_c_api.h>
#include <test_config.hpp>

int main()
{
   std::string file = TESTS_PATH "/example.config";

   set_up_(file.c_str(), 0);
   return 0;
}