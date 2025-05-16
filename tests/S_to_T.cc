
#include <honeycomb2/honeycomb2.hpp>
#include "discretization.hpp"
#include "honeycomb_120_25_grid_points.hpp"
#include "original_model_functions.hpp"
#include "solution.hpp"
#include "utilities.hpp"

int main()
{

   const size_t n         = 8;
   const double rmin      = 0.01;
   Honeycomb::Grid2D grid = Honeycomb::generate_compliant_Grid2D(n, {rmin, 0.15, 0.65, 1}, {13, 10, 7});

   Honeycomb::Discretization discr(grid);

   Honeycomb::InputModel model;

   model.SetModel(Honeycomb::InputModel::T_UP, orignal_models::Tu_test);
   model.SetModel(Honeycomb::InputModel::DT_UP, orignal_models::DTu_test);

   model.SetModel(Honeycomb::InputModel::T_DN, orignal_models::Td_test);
   model.SetModel(Honeycomb::InputModel::DT_DN, orignal_models::DTd_test);

   model.SetModel(Honeycomb::InputModel::T_ST, orignal_models::Ts_test);
   model.SetModel(Honeycomb::InputModel::DT_ST, orignal_models::DTs_test);

   model.SetModel(Honeycomb::InputModel::T_P_GL, orignal_models::TFp_test);
   model.SetModel(Honeycomb::InputModel::T_M_GL, orignal_models::TFm_test);

   Honeycomb::Solution sol0(&discr, model, 3);
   sol0.RotateToPhysicalBasis();

   Honeycomb::OutputModel o_model(sol0);

   double m = 0;

   for (long int i = 0; i < grid.c_size_li; i++) {
      Honeycomb::RnC::Triplet x123 = grid._x123[i];

      m = std::max(m, std::fabs(orignal_models::TFp_test(x123(0), x123(1), x123(2)) - o_model.T[0][i]));
      m = std::max(m, std::fabs(orignal_models::Td_test(x123(0), x123(1), x123(2)) - o_model.T[1][i]));
      m = std::max(m, std::fabs(orignal_models::Tu_test(x123(0), x123(1), x123(2)) - o_model.T[2][i]));
      m = std::max(m, std::fabs(orignal_models::Ts_test(x123(0), x123(1), x123(2)) - o_model.T[3][i]));
      m = std::max(
          m, std::fabs(orignal_models::model_zero_function(x123(0), x123(1), x123(2)) - o_model.T[4][i]));
      m = std::max(
          m, std::fabs(orignal_models::model_zero_function(x123(0), x123(1), x123(2)) - o_model.T[5][i]));
      m = std::max(
          m, std::fabs(orignal_models::model_zero_function(x123(0), x123(1), x123(2)) - o_model.T[6][i]));

      m = std::max(m, std::fabs(orignal_models::TFm_test(x123(0), x123(1), x123(2)) - o_model.DT[0][i]));
      m = std::max(m, std::fabs(orignal_models::DTd_test(x123(0), x123(1), x123(2)) - o_model.DT[1][i]));
      m = std::max(m, std::fabs(orignal_models::DTu_test(x123(0), x123(1), x123(2)) - o_model.DT[2][i]));
      m = std::max(m, std::fabs(orignal_models::DTs_test(x123(0), x123(1), x123(2)) - o_model.DT[3][i]));
      m = std::max(
          m, std::fabs(orignal_models::model_zero_function(x123(0), x123(1), x123(2)) - o_model.DT[4][i]));
      m = std::max(
          m, std::fabs(orignal_models::model_zero_function(x123(0), x123(1), x123(2)) - o_model.DT[5][i]));
      m = std::max(
          m, std::fabs(orignal_models::model_zero_function(x123(0), x123(1), x123(2)) - o_model.DT[6][i]));
   }

   Honeycomb::logger(Honeycomb::Logger::WARNING, std::format("Max difference: {:.10e}", m));

   m = 0;
   for (const Honeycomb::RnC::Triplet &x123 : honeycomb_points) {
      double exact  = 0;
      double approx = 0;

      exact  = orignal_models::TFp_test(x123(0), x123(1), x123(2));
      approx = o_model.GetDistribution(Honeycomb::OutputModel::T_P_GL, x123);
      m      = std::max(m, std::fabs(exact - approx));

      exact  = orignal_models::Td_test(x123(0), x123(1), x123(2));
      approx = o_model.GetDistribution(Honeycomb::OutputModel::T_DN, x123);
      m      = std::max(m, std::fabs(exact - approx));

      exact  = orignal_models::Tu_test(x123(0), x123(1), x123(2));
      approx = o_model.GetDistribution(Honeycomb::OutputModel::T_UP, x123);
      m      = std::max(m, std::fabs(exact - approx));

      exact  = orignal_models::Ts_test(x123(0), x123(1), x123(2));
      approx = o_model.GetDistribution(Honeycomb::OutputModel::T_ST, x123);
      m      = std::max(m, std::fabs(exact - approx));

      exact  = orignal_models::model_zero_function(x123(0), x123(1), x123(2));
      approx = o_model.GetDistribution(Honeycomb::OutputModel::T_CH, x123);
      m      = std::max(m, std::fabs(exact - approx));

      exact  = orignal_models::model_zero_function(x123(0), x123(1), x123(2));
      approx = o_model.GetDistribution(Honeycomb::OutputModel::T_BM, x123);
      m      = std::max(m, std::fabs(exact - approx));

      exact  = orignal_models::model_zero_function(x123(0), x123(1), x123(2));
      approx = o_model.GetDistribution(Honeycomb::OutputModel::T_TP, x123);
      m      = std::max(m, std::fabs(exact - approx));

      // -----------------------------------------------------------------------

      exact  = orignal_models::TFm_test(x123(0), x123(1), x123(2));
      approx = o_model.GetDistribution(Honeycomb::OutputModel::T_M_GL, x123);
      m      = std::max(m, std::fabs(exact - approx));

      exact  = orignal_models::DTd_test(x123(0), x123(1), x123(2));
      approx = o_model.GetDistribution(Honeycomb::OutputModel::DT_DN, x123);
      m      = std::max(m, std::fabs(exact - approx));

      exact  = orignal_models::DTu_test(x123(0), x123(1), x123(2));
      approx = o_model.GetDistribution(Honeycomb::OutputModel::DT_UP, x123);
      m      = std::max(m, std::fabs(exact - approx));

      exact  = orignal_models::DTs_test(x123(0), x123(1), x123(2));
      approx = o_model.GetDistribution(Honeycomb::OutputModel::DT_ST, x123);
      m      = std::max(m, std::fabs(exact - approx));

      exact  = orignal_models::model_zero_function(x123(0), x123(1), x123(2));
      approx = o_model.GetDistribution(Honeycomb::OutputModel::DT_CH, x123);
      m      = std::max(m, std::fabs(exact - approx));

      exact  = orignal_models::model_zero_function(x123(0), x123(1), x123(2));
      approx = o_model.GetDistribution(Honeycomb::OutputModel::DT_BM, x123);
      m      = std::max(m, std::fabs(exact - approx));

      exact  = orignal_models::model_zero_function(x123(0), x123(1), x123(2));
      approx = o_model.GetDistribution(Honeycomb::OutputModel::DT_TP, x123);
      m      = std::max(m, std::fabs(exact - approx));
   }

   Honeycomb::logger(Honeycomb::Logger::WARNING, std::format("Max difference: {:.10e}", m));
}