
#include <honeycomb2/honeycomb2.hpp>
#include "honeycomb_120_25_grid_points.hpp"
#include "original_model_functions.hpp"

int main()
{

   const size_t n    = 8;
   const double rmin = 0.01;
   Honeycomb::Grid2D grid
       = Honeycomb::generate_compliant_Grid2D(n, {rmin, 0.15, 0.65, 1}, {13, 10, 7});

   Honeycomb::Discretization discr(grid);

   Honeycomb::InputModel model = Honeycomb::PreImplementedModels::GetModel("pim_original");

   Honeycomb::Solution sol0(&discr, model, 3);
   sol0.RotateToPhysicalBasis();

   Honeycomb::OutputModel o_model(sol0);

   double m = 0;

   for (long int i = 0; i < grid.c_size_li; i++) {
      Honeycomb::RnC::Triplet x123 = grid._x123[i];

      m = std::fabs(orignal_models::TFp_test(x123(0), x123(1), x123(2)) - o_model.T[0][i]);
      if (!Honeycomb::is_near(m, 0.0, 1.0e-14)) {
         std::cout << m << std::endl;
         return -1;
      }
      m = std::fabs(orignal_models::Td_test(x123(0), x123(1), x123(2)) - o_model.T[1][i]);
      if (!Honeycomb::is_near(m, 0.0, 1.0e-14)) {
         std::cout << m << std::endl;
         return -1;
      }
      m = std::fabs(orignal_models::Tu_test(x123(0), x123(1), x123(2)) - o_model.T[2][i]);
      if (!Honeycomb::is_near(m, 0.0, 1.0e-14)) {
         std::cout << m << std::endl;
         return -1;
      }
      m = std::fabs(orignal_models::Ts_test(x123(0), x123(1), x123(2)) - o_model.T[3][i]);
      if (!Honeycomb::is_near(m, 0.0, 1.0e-14)) {
         std::cout << m << std::endl;
         return -1;
      }
      m = std::fabs(orignal_models::model_zero_function(x123(0), x123(1), x123(2))
                    - o_model.T[4][i]);
      if (!Honeycomb::is_near(m, 0.0, 1.0e-14)) {
         std::cout << m << std::endl;
         return -1;
      }
      m = std::fabs(orignal_models::model_zero_function(x123(0), x123(1), x123(2))
                    - o_model.T[5][i]);
      if (!Honeycomb::is_near(m, 0.0, 1.0e-14)) {
         std::cout << m << std::endl;
         return -1;
      }
      m = std::fabs(orignal_models::model_zero_function(x123(0), x123(1), x123(2))
                    - o_model.T[6][i]);
      if (!Honeycomb::is_near(m, 0.0, 1.0e-14)) {
         std::cout << m << std::endl;
         return -1;
      }

      m = std::fabs(orignal_models::TFm_test(x123(0), x123(1), x123(2)) - o_model.DT[0][i]);
      if (!Honeycomb::is_near(m, 0.0, 1.0e-14)) {
         std::cout << m << std::endl;
         return -1;
      }
      m = std::fabs(orignal_models::DTd_test(x123(0), x123(1), x123(2)) - o_model.DT[1][i]);
      if (!Honeycomb::is_near(m, 0.0, 1.0e-14)) {
         std::cout << m << std::endl;
         return -1;
      }
      m = std::fabs(orignal_models::DTu_test(x123(0), x123(1), x123(2)) - o_model.DT[2][i]);
      if (!Honeycomb::is_near(m, 0.0, 1.0e-14)) {
         std::cout << m << std::endl;
         return -1;
      }
      m = std::fabs(orignal_models::DTs_test(x123(0), x123(1), x123(2)) - o_model.DT[3][i]);
      if (!Honeycomb::is_near(m, 0.0, 1.0e-14)) {
         std::cout << m << std::endl;
         return -1;
      }
      m = std::fabs(orignal_models::model_zero_function(x123(0), x123(1), x123(2))
                    - o_model.DT[4][i]);
      if (!Honeycomb::is_near(m, 0.0, 1.0e-14)) {
         std::cout << m << std::endl;
         return -1;
      }
      m = std::fabs(orignal_models::model_zero_function(x123(0), x123(1), x123(2))
                    - o_model.DT[5][i]);
      if (!Honeycomb::is_near(m, 0.0, 1.0e-14)) {
         std::cout << m << std::endl;
         return -1;
      }
      m = std::fabs(orignal_models::model_zero_function(x123(0), x123(1), x123(2))
                    - o_model.DT[6][i]);
      if (!Honeycomb::is_near(m, 0.0, 1.0e-14)) {
         std::cout << m << std::endl;
         return -1;
      }
   }
   return 0;
}