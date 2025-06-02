#include <honeycomb2/honeycomb2.hpp>
#include <honeycomb2/honeycomb2_c_api.h>
#include <test_config.hpp>
#include "honeycomb_120_25_grid_points.hpp"
#include "original_model_functions.hpp"

double Tu_test_loc(double *x1, double *x2, double *x3)
{
   return orignal_models::Tu_test(*x1, *x2, *x3);
};
double DTu_test_loc(double *x1, double *x2, double *x3)
{
   return orignal_models::DTu_test(*x1, *x2, *x3);
};
double Td_test_loc(double *x1, double *x2, double *x3)
{
   return orignal_models::Td_test(*x1, *x2, *x3);
};
double DTd_test_loc(double *x1, double *x2, double *x3)
{
   return orignal_models::DTd_test(*x1, *x2, *x3);
};
double Ts_test_loc(double *x1, double *x2, double *x3)
{
   return orignal_models::Ts_test(*x1, *x2, *x3);
};
double DTs_test_loc(double *x1, double *x2, double *x3)
{
   return orignal_models::DTs_test(*x1, *x2, *x3);
};
double TFp_test_loc(double *x1, double *x2, double *x3)
{
   return orignal_models::TFp_test(*x1, *x2, *x3);
};
double TFm_test_loc(double *x1, double *x2, double *x3)
{
   return orignal_models::TFm_test(*x1, *x2, *x3);
};

int main()
{
   std::string file = TESTS_PATH "/example.config";

   hc2_fi_set_up_(file.c_str(), file.size());

   int i = 0;

   i = static_cast<int>(Honeycomb::InputModel::T_UP);
   hc2_fi_set_model_(&i, Tu_test_loc);
   i = static_cast<int>(Honeycomb::InputModel::DT_UP);
   hc2_fi_set_model_(&i, DTu_test_loc);
   i = static_cast<int>(Honeycomb::InputModel::T_DN);
   hc2_fi_set_model_(&i, Td_test_loc);
   i = static_cast<int>(Honeycomb::InputModel::DT_DN);
   hc2_fi_set_model_(&i, DTd_test_loc);
   i = static_cast<int>(Honeycomb::InputModel::T_ST);
   hc2_fi_set_model_(&i, Ts_test_loc);
   i = static_cast<int>(Honeycomb::InputModel::DT_ST);
   hc2_fi_set_model_(&i, DTs_test_loc);
   i = static_cast<int>(Honeycomb::InputModel::T_P_GL);
   hc2_fi_set_model_(&i, TFp_test_loc);
   i = static_cast<int>(Honeycomb::InputModel::T_M_GL);
   hc2_fi_set_model_(&i, TFm_test_loc);

   hc2_fi_evolve_();
   hc2_fi_unload_();

   double Q2fin = 4;
   double Q20   = 1;

   std::FILE *fp   = std::fopen("T_and_DT_final_scale.dat", "w");
   std::FILE *fp_1 = std::fopen("T_and_DT_initial_scale.dat", "w");
   for (Honeycomb::RnC::Triplet &x123 : honeycomb_points) {

      // Column 1, 2, and 3 are x1, x2, x3
      std::fprintf(fp, "%.16e\t%.16e\t%.16e\t", x123(0), x123(1), x123(2));
      std::fprintf(fp_1, "%.16e\t%.16e\t%.16e\t", x123(0), x123(1), x123(2));

      // Save distributions. Order is:
      // T_{3F}^+, T_{3F}^-, T_d, DT_d, T_u, DT_u, ...
      // ^^^^^^^^^^^^^^^^^^
      // These are in Eq. 23 of 2405.01162

      double x1 = x123(0);
      double x2 = x123(1);
      double x3 = x123(2);

      //------------------------------------------
      // This is final scale model
      i = static_cast<int>(Honeycomb::OutputModel::T_P_GL);
      std::fprintf(fp, "%.16e\t", hc2_fi_get_model_(&i, &Q2fin, &x1, &x2, &x3));
      i = static_cast<int>(Honeycomb::OutputModel::T_M_GL);
      std::fprintf(fp, "%.16e\t", hc2_fi_get_model_(&i, &Q2fin, &x1, &x2, &x3));
      i = static_cast<int>(Honeycomb::OutputModel::T_DN);
      std::fprintf(fp, "%.16e\t", hc2_fi_get_model_(&i, &Q2fin, &x1, &x2, &x3));
      i = static_cast<int>(Honeycomb::OutputModel::DT_DN);
      std::fprintf(fp, "%.16e\t", hc2_fi_get_model_(&i, &Q2fin, &x1, &x2, &x3));
      i = static_cast<int>(Honeycomb::OutputModel::T_UP);
      std::fprintf(fp, "%.16e\t", hc2_fi_get_model_(&i, &Q2fin, &x1, &x2, &x3));
      i = static_cast<int>(Honeycomb::OutputModel::DT_UP);
      std::fprintf(fp, "%.16e\t", hc2_fi_get_model_(&i, &Q2fin, &x1, &x2, &x3));
      i = static_cast<int>(Honeycomb::OutputModel::T_ST);
      std::fprintf(fp, "%.16e\t", hc2_fi_get_model_(&i, &Q2fin, &x1, &x2, &x3));
      i = static_cast<int>(Honeycomb::OutputModel::DT_ST);
      std::fprintf(fp, "%.16e\t", hc2_fi_get_model_(&i, &Q2fin, &x1, &x2, &x3));
      i = static_cast<int>(Honeycomb::OutputModel::T_CH);
      std::fprintf(fp, "%.16e\t", hc2_fi_get_model_(&i, &Q2fin, &x1, &x2, &x3));
      i = static_cast<int>(Honeycomb::OutputModel::DT_CH);
      std::fprintf(fp, "%.16e\t", hc2_fi_get_model_(&i, &Q2fin, &x1, &x2, &x3));
      i = static_cast<int>(Honeycomb::OutputModel::T_BM);
      std::fprintf(fp, "%.16e\t", hc2_fi_get_model_(&i, &Q2fin, &x1, &x2, &x3));
      i = static_cast<int>(Honeycomb::OutputModel::DT_BM);
      std::fprintf(fp, "%.16e\t", hc2_fi_get_model_(&i, &Q2fin, &x1, &x2, &x3));

      // Omit Top quark, essentially irrelevant for us.
      std::fprintf(fp, "\n");

      //------------------------------------------
      // This is initial scale model
      i = static_cast<int>(Honeycomb::OutputModel::T_P_GL);
      std::fprintf(fp_1, "%.16e\t", hc2_fi_get_model_(&i, &Q20, &x1, &x2, &x3));
      i = static_cast<int>(Honeycomb::OutputModel::T_M_GL);
      std::fprintf(fp_1, "%.16e\t", hc2_fi_get_model_(&i, &Q20, &x1, &x2, &x3));
      i = static_cast<int>(Honeycomb::OutputModel::T_DN);
      std::fprintf(fp_1, "%.16e\t", hc2_fi_get_model_(&i, &Q20, &x1, &x2, &x3));
      i = static_cast<int>(Honeycomb::OutputModel::DT_DN);
      std::fprintf(fp_1, "%.16e\t", hc2_fi_get_model_(&i, &Q20, &x1, &x2, &x3));
      i = static_cast<int>(Honeycomb::OutputModel::T_UP);
      std::fprintf(fp_1, "%.16e\t", hc2_fi_get_model_(&i, &Q20, &x1, &x2, &x3));
      i = static_cast<int>(Honeycomb::OutputModel::DT_UP);
      std::fprintf(fp_1, "%.16e\t", hc2_fi_get_model_(&i, &Q20, &x1, &x2, &x3));
      i = static_cast<int>(Honeycomb::OutputModel::T_ST);
      std::fprintf(fp_1, "%.16e\t", hc2_fi_get_model_(&i, &Q20, &x1, &x2, &x3));
      i = static_cast<int>(Honeycomb::OutputModel::DT_ST);
      std::fprintf(fp_1, "%.16e\t", hc2_fi_get_model_(&i, &Q20, &x1, &x2, &x3));
      i = static_cast<int>(Honeycomb::OutputModel::T_CH);
      std::fprintf(fp_1, "%.16e\t", hc2_fi_get_model_(&i, &Q20, &x1, &x2, &x3));
      i = static_cast<int>(Honeycomb::OutputModel::DT_CH);
      std::fprintf(fp_1, "%.16e\t", hc2_fi_get_model_(&i, &Q20, &x1, &x2, &x3));
      i = static_cast<int>(Honeycomb::OutputModel::T_BM);
      std::fprintf(fp_1, "%.16e\t", hc2_fi_get_model_(&i, &Q20, &x1, &x2, &x3));
      i = static_cast<int>(Honeycomb::OutputModel::DT_BM);
      std::fprintf(fp_1, "%.16e\t", hc2_fi_get_model_(&i, &Q20, &x1, &x2, &x3));

      // Omit Top quark, essentially irrelevant for us.
      std::fprintf(fp_1, "\n");
   }

   return 0;
}