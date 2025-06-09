//
// Author(s): Simone Rodini   <rodini.simone.luigi@gmail.com>
//            Arianna Vercesi <arianna.vercesi01@universitadipavia.it>
//

#include <honeycomb2/honeycomb2.hpp>
#include "test_rd_LWFA_asymptotic.hpp"
#include "test_rd_LWFA_fitted.hpp"

int main()
{
   if (!Honeycomb::PreImplementedModels::QueryModelAvailability("pim_original")) return 1;
   if (!Honeycomb::PreImplementedModels::QueryModelAvailability("pim_LFWA_asymptotic")) return 1;
   if (!Honeycomb::PreImplementedModels::QueryModelAvailability("pim_LFWA_fitted")) return 1;

   Honeycomb::InputModel model = Honeycomb::PreImplementedModels::GetModel("pim_LFWA_asymptotic");
   for (const std::vector<double> &entry : T_DT_LFWA_asymptotic) {
      double T_dn = model.T[0](entry[0], entry[1], entry[2]);
      double T_up = model.T[1](entry[0], entry[1], entry[2]);

      double DT_dn = model.DT[0](entry[0], entry[1], entry[2]);
      double DT_up = model.DT[1](entry[0], entry[1], entry[2]);

      if (!Honeycomb::is_near(T_dn, entry[3], 1.0e-12)) return 1;
      if (!Honeycomb::is_near(DT_dn, entry[4], 1.0e-12)) return 1;
      if (!Honeycomb::is_near(T_up, entry[5], 1.0e-12)) return 1;
      if (!Honeycomb::is_near(DT_up, entry[6], 1.0e-12)) return 1;
   }

   model = Honeycomb::PreImplementedModels::GetModel("pim_LFWA_fitted");
   for (const std::vector<double> &entry : T_DT_LFWA_fitted) {
      double T_dn = model.T[0](entry[0], entry[1], entry[2]);
      double T_up = model.T[1](entry[0], entry[1], entry[2]);

      double DT_dn = model.DT[0](entry[0], entry[1], entry[2]);
      double DT_up = model.DT[1](entry[0], entry[1], entry[2]);

      if (!Honeycomb::is_near(T_dn, entry[3], 1.0e-12)) return 1;
      if (!Honeycomb::is_near(DT_dn, entry[4], 1.0e-12)) return 1;
      if (!Honeycomb::is_near(T_up, entry[5], 1.0e-12)) return 1;
      if (!Honeycomb::is_near(DT_up, entry[6], 1.0e-12)) return 1;
   }

   return 0;
}