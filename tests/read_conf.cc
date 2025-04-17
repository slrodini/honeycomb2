#include "honeycomb2/config_parser.hpp"
#include <honeycomb2/honeycomb2.hpp>

int main()
{
   //
   std::string content = "key: \"val kappa\"\nFoo: 4.2";
   Honeycomb::ConfigParser parser(content);

   double q          = parser.GetValue<double>("Foo");
   std::string check = parser.GetValue("key");
   std::cout << "|" << check << "|" << std::endl;

   std::cout << q << std::endl;
   return 0;
}