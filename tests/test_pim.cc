#include "solution.hpp"
#include <honeycomb2/honeycomb2.hpp>

int main()
{
   if (!Honeycomb::PreImplementedModels::QueryModelAvailability("pim_original")) return 1;

   if (!Honeycomb::PreImplementedModels::QueryModelAvailability("pim_LFWA_asymptotic")) return 1;

   if (!Honeycomb::PreImplementedModels::QueryModelAvailability("pim_LFWA_fitted")) return 1;

   return 0;
}