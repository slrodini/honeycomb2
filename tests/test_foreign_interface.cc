#include <honeycomb2/honeycomb2.hpp>
#include <honeycomb2/honeycomb2_c_api.h>
#include <test_config.hpp>
#include "honeycomb_120_25_grid_points.hpp"
#include "test_rd_T_and_DT_final_scale.hpp"
#include "test_rd_T_and_DT_initial_scale.hpp"

int main()
{
   std::string file = TESTS_PATH "/example.config";

   hc2_fi_set_up_(file.c_str(), file.size());

   hc2_fi_evolve_();
   hc2_fi_unload_();

   double Q2fin = 10000;
   double Q20   = 1;

   std::FILE *fp   = std::fopen("T_and_DT_final_scale.dat", "w");
   std::FILE *fp_1 = std::fopen("T_and_DT_initial_scale.dat", "w");

   int i;
   size_t index = 0;
   for (Honeycomb::RnC::Triplet &x123 : honeycomb_points) {

      // Column 1, 2, and 3 are x1, x2, x3
      std::fprintf(fp, "%.16e\t%.16e\t%.16e\t", x123(0), x123(1), x123(2));
      std::fprintf(fp_1, "%.16e\t%.16e\t%.16e\t", x123(0), x123(1), x123(2));

      // Save distributions. Order is:
      // T_{3F}^+, T_{3F}^-, T_d, DT_d, T_u, DT_u, ...
      // ^^^^^^^^^^^^^^^^^^
      // These are in Eq. 23 of 2405.01162

      double x1      = x123(0);
      double x2      = x123(1);
      double x3      = x123(2);
      double curr_pt = 0;
      int ret_val    = 0;

      //------------------------------------------
      // This is final scale model
      i       = static_cast<int>(Honeycomb::OutputModel::T_P_GL);
      curr_pt = hc2_fi_get_model_(&i, &Q2fin, &x1, &x2, &x3);
      if (!Honeycomb::is_near(curr_pt, T_DT_regression_final[index][3], 1.0e-6)) ret_val += 1;

      i       = static_cast<int>(Honeycomb::OutputModel::T_M_GL);
      curr_pt = hc2_fi_get_model_(&i, &Q2fin, &x1, &x2, &x3);
      if (!Honeycomb::is_near(curr_pt, T_DT_regression_final[index][4], 1.0e-6)) ret_val += 1;

      i       = static_cast<int>(Honeycomb::OutputModel::T_DN);
      curr_pt = hc2_fi_get_model_(&i, &Q2fin, &x1, &x2, &x3);
      if (!Honeycomb::is_near(curr_pt, T_DT_regression_final[index][5], 1.0e-6)) ret_val += 1;

      i       = static_cast<int>(Honeycomb::OutputModel::DT_DN);
      curr_pt = hc2_fi_get_model_(&i, &Q2fin, &x1, &x2, &x3);
      if (!Honeycomb::is_near(curr_pt, T_DT_regression_final[index][6], 1.0e-6)) ret_val += 1;

      i       = static_cast<int>(Honeycomb::OutputModel::T_UP);
      curr_pt = hc2_fi_get_model_(&i, &Q2fin, &x1, &x2, &x3);
      if (!Honeycomb::is_near(curr_pt, T_DT_regression_final[index][7], 1.0e-6)) ret_val += 1;

      i       = static_cast<int>(Honeycomb::OutputModel::DT_UP);
      curr_pt = hc2_fi_get_model_(&i, &Q2fin, &x1, &x2, &x3);
      if (!Honeycomb::is_near(curr_pt, T_DT_regression_final[index][8], 1.0e-6)) ret_val += 1;

      i       = static_cast<int>(Honeycomb::OutputModel::T_ST);
      curr_pt = hc2_fi_get_model_(&i, &Q2fin, &x1, &x2, &x3);
      if (!Honeycomb::is_near(curr_pt, T_DT_regression_final[index][9], 1.0e-6)) ret_val += 1;

      i       = static_cast<int>(Honeycomb::OutputModel::DT_ST);
      curr_pt = hc2_fi_get_model_(&i, &Q2fin, &x1, &x2, &x3);
      if (!Honeycomb::is_near(curr_pt, T_DT_regression_final[index][10], 1.0e-6)) ret_val += 1;

      i       = static_cast<int>(Honeycomb::OutputModel::T_CH);
      curr_pt = hc2_fi_get_model_(&i, &Q2fin, &x1, &x2, &x3);
      if (!Honeycomb::is_near(curr_pt, T_DT_regression_final[index][11], 1.0e-6)) ret_val += 1;

      i       = static_cast<int>(Honeycomb::OutputModel::DT_CH);
      curr_pt = hc2_fi_get_model_(&i, &Q2fin, &x1, &x2, &x3);
      if (!Honeycomb::is_near(curr_pt, T_DT_regression_final[index][12], 1.0e-6)) ret_val += 1;

      i       = static_cast<int>(Honeycomb::OutputModel::T_BM);
      curr_pt = hc2_fi_get_model_(&i, &Q2fin, &x1, &x2, &x3);
      if (!Honeycomb::is_near(curr_pt, T_DT_regression_final[index][13], 1.0e-6)) ret_val += 1;

      i       = static_cast<int>(Honeycomb::OutputModel::DT_BM);
      curr_pt = hc2_fi_get_model_(&i, &Q2fin, &x1, &x2, &x3);
      if (!Honeycomb::is_near(curr_pt, T_DT_regression_final[index][14], 1.0e-6)) ret_val += 1;

      //------------------------------------------
      // This is initial scale model
      i       = static_cast<int>(Honeycomb::OutputModel::T_P_GL);
      curr_pt = hc2_fi_get_model_(&i, &Q20, &x1, &x2, &x3);
      if (!Honeycomb::is_near(curr_pt, T_DT_regression_initial[index][3], 1.0e-6)) ret_val += 1;

      i       = static_cast<int>(Honeycomb::OutputModel::T_M_GL);
      curr_pt = hc2_fi_get_model_(&i, &Q20, &x1, &x2, &x3);
      if (!Honeycomb::is_near(curr_pt, T_DT_regression_initial[index][4], 1.0e-6)) ret_val += 1;

      i       = static_cast<int>(Honeycomb::OutputModel::T_DN);
      curr_pt = hc2_fi_get_model_(&i, &Q20, &x1, &x2, &x3);
      if (!Honeycomb::is_near(curr_pt, T_DT_regression_initial[index][5], 1.0e-6)) ret_val += 1;

      i       = static_cast<int>(Honeycomb::OutputModel::DT_DN);
      curr_pt = hc2_fi_get_model_(&i, &Q20, &x1, &x2, &x3);
      if (!Honeycomb::is_near(curr_pt, T_DT_regression_initial[index][6], 1.0e-6)) ret_val += 1;

      i       = static_cast<int>(Honeycomb::OutputModel::T_UP);
      curr_pt = hc2_fi_get_model_(&i, &Q20, &x1, &x2, &x3);
      if (!Honeycomb::is_near(curr_pt, T_DT_regression_initial[index][7], 1.0e-6)) ret_val += 1;

      i       = static_cast<int>(Honeycomb::OutputModel::DT_UP);
      curr_pt = hc2_fi_get_model_(&i, &Q20, &x1, &x2, &x3);
      if (!Honeycomb::is_near(curr_pt, T_DT_regression_initial[index][8], 1.0e-6)) ret_val += 1;

      i       = static_cast<int>(Honeycomb::OutputModel::T_ST);
      curr_pt = hc2_fi_get_model_(&i, &Q20, &x1, &x2, &x3);
      if (!Honeycomb::is_near(curr_pt, T_DT_regression_initial[index][9], 1.0e-6)) ret_val += 1;

      i       = static_cast<int>(Honeycomb::OutputModel::DT_ST);
      curr_pt = hc2_fi_get_model_(&i, &Q20, &x1, &x2, &x3);
      if (!Honeycomb::is_near(curr_pt, T_DT_regression_initial[index][10], 1.0e-6)) ret_val += 1;

      i       = static_cast<int>(Honeycomb::OutputModel::T_CH);
      curr_pt = hc2_fi_get_model_(&i, &Q20, &x1, &x2, &x3);
      if (!Honeycomb::is_near(curr_pt, T_DT_regression_initial[index][11], 1.0e-6)) ret_val += 1;

      i       = static_cast<int>(Honeycomb::OutputModel::DT_CH);
      curr_pt = hc2_fi_get_model_(&i, &Q20, &x1, &x2, &x3);
      if (!Honeycomb::is_near(curr_pt, T_DT_regression_initial[index][12], 1.0e-6)) ret_val += 1;

      i       = static_cast<int>(Honeycomb::OutputModel::T_BM);
      curr_pt = hc2_fi_get_model_(&i, &Q20, &x1, &x2, &x3);
      if (!Honeycomb::is_near(curr_pt, T_DT_regression_initial[index][13], 1.0e-6)) ret_val += 1;

      i       = static_cast<int>(Honeycomb::OutputModel::DT_BM);
      curr_pt = hc2_fi_get_model_(&i, &Q20, &x1, &x2, &x3);
      if (!Honeycomb::is_near(curr_pt, T_DT_regression_initial[index][14], 1.0e-6)) ret_val += 1;

      if (ret_val != 0) return ret_val;

      index++;
   }

   return 0;
}