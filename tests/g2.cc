#include <honeycomb2/honeycomb2.hpp>
#include "honeycomb2/obs_weights.hpp"
#include "honeycomb2/solution.hpp"
#include "honeycomb_120_25_grid_points.hpp"
#include <cereal/archives/json.hpp>

#include "original_model_functions.hpp"

inline void print_grid_dim(const Honeycomb::Grid2D &grid_1)
{
   Honeycomb::logger(Honeycomb::Logger::INFO,
                     std::format("Total grid size: {:d}x{:d}", grid_1.size, grid_1.size));
   Honeycomb::logger(Honeycomb::Logger::INFO,
                     std::format("Total x123 size: {:d}x{:d}", grid_1.c_size, grid_1.c_size));
}

void compute_one_grid(const Honeycomb::Grid2D &grid, const std::string file_name);

int main()
{
   //
   const double rmin = 0.01;

   // Honeycomb::Grid2D grid_1 = Honeycomb::generate_compliant_Grid2D(2, {rmin, 0.4, 1}, {6, 6});
   // print_grid_dim(grid_1);

   // Honeycomb::Grid2D grid_2 = Honeycomb::generate_compliant_Grid2D(4, {rmin, 0.4, 1}, {10, 10});
   // print_grid_dim(grid_2);

   // Honeycomb::Grid2D grid_3 = Honeycomb::generate_compliant_Grid2D(4, {rmin, 0.1, 0.4, 1}, {12,
   // 9, 7}); print_grid_dim(grid_3);

   // Honeycomb::Grid2D grid_4 = Honeycomb::generate_compliant_Grid2D(8, {rmin, 0.1, 0.4, 1}, {12,
   // 9, 7}); print_grid_dim(grid_4);

   // Honeycomb::Grid2D grid_5 = Honeycomb::generate_compliant_Grid2D(20, {rmin, 0.1, 0.4, 1}, {12,
   // 9, 7}); print_grid_dim(grid_5);

   Honeycomb::Grid2D grid_6
       = Honeycomb::generate_compliant_Grid2D(8, {rmin, 0.1, 0.4, 1}, {12, 9, 7});
   print_grid_dim(grid_6);

   // compute_one_grid(grid_1, "g2_grid_1.dat");
   // compute_one_grid(grid_2, "g2_grid_2.dat");
   // compute_one_grid(grid_3, "g2_grid_3.dat");
   // compute_one_grid(grid_4, "g2_grid_4.dat");
   // compute_one_grid(grid_5, "g2_grid_5.dat");
   compute_one_grid(grid_6, "g2_grid_6.dat");

   return 0;
}

void compute_one_grid(const Honeycomb::Grid2D &grid, const std::string file_name)
{
   Honeycomb::logger(Honeycomb::Logger::INFO,
                     std::format("================== Begin: {:s} =================", file_name));

   double fnc_elapsed           = 0;
   Honeycomb::timer::mark begin = Honeycomb::timer::now();
   Honeycomb::timer::mark end   = Honeycomb::timer::now();

   const double Nc = 3; // NC = 1 for tests

   Honeycomb::logger(Honeycomb::Logger::INFO, std::format("Kernel Discretization"));
   begin = Honeycomb::timer::now();
   Honeycomb::Kernels kers(grid, Nc);

   end = Honeycomb::timer::now();
   Honeycomb::logger(
       Honeycomb::Logger::INFO,
       std::format("  Elapsed: {:.4e} (ms)", Honeycomb::timer::elapsed_ms(end, begin)));

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

   auto as = [](double t) -> double {
      return 1.0 / (11.0 * (t + 3.0));
   };
   const double pref = -1.0;

   const double Q02 = 1.0;
   const double Qf2 = 1.0e+2;

   const double t0 = log(Q02);

   auto [inter_scale, sol1]
       = get_initial_solution(Q02, Qf2, {0, 0, 0, 1.6129, 17.4724, 1.0e+6}, &discr, model);

   const double dt = 0.01;

   Honeycomb::runge_kutta::GenericRungeKutta<Honeycomb::Kernels, Honeycomb::Solution, 13> evolver(
       kers, sol1, Honeycomb::runge_kutta::DOPRI8, as, pref, t0, dt);

   std::vector<size_t> steps(inter_scale.size(), 20);
   fnc_elapsed = 0;
   begin       = Honeycomb::timer::now();
   evolver(inter_scale, steps);
   end         = Honeycomb::timer::now();
   fnc_elapsed = Honeycomb::timer::elapsed_ms(end, begin);

   Honeycomb::logger(Honeycomb::Logger::INFO, std::format("  Elapsed: {:.4e} (ms)", fnc_elapsed));

   auto sol_fin = evolver.GetSolution();

   sol_fin.RotateToPhysicalBasis();

   {
      Eigen::VectorXd Splus_summed = Eigen::VectorXd::Zero(sol_fin._distr_p[0].size());
      for (size_t j = 1; j < sol_fin._distr_p.size(); j++) {
         double pref = 1.0 / 9.0;
         if (j % 2 == 0) pref *= 4; // Up-type quarks are j=2, 4, 6
         Splus_summed += pref * sol_fin._distr_p[j];
      }

      std::FILE *fp_3
          = std::fopen(std::string("CE_evo_check_evo_basis_" + file_name + ".dat").c_str(), "w");

      for (Honeycomb::RnC::Triplet &x123 : honeycomb_points) {
         Honeycomb::RnC::Pair rhophi = Honeycomb::RnC::from_x123_to_rhophi(x123);

         std::fprintf(fp_3, "%.16e\t%.16e\t%.16e\t", x123[0], x123[1], x123[2]);
         std::fprintf(fp_3, "%.16e\t", discr.interpolate_as_weights_v3(rhophi, Splus_summed));
         std::fprintf(fp_3, "\n");
      }
      std::fclose(fp_3);
   }

   std::vector<Honeycomb::G2Weights> weights;
   std::vector<double> xBjs = {0.1,  0.15, 0.2,  0.25, 0.3,  0.35, 0.4,  0.45, 0.5,
                               0.55, 0.6,  0.65, 0.7,  0.75, 0.8,  0.85, 0.9,  0.95};

   double total_time = 0;
   for (const double &xBj : xBjs) {
      begin = Honeycomb::timer::now();
      weights.emplace_back(Honeycomb::G2Weights(xBj, grid));
      end         = Honeycomb::timer::now();
      fnc_elapsed = Honeycomb::timer::elapsed_ms(end, begin);

      Honeycomb::logger(Honeycomb::Logger::INFO,
                        std::format("  xBj: {:2f}; Elapsed: {:.4e} (ms)", xBj, fnc_elapsed));
      total_time += fnc_elapsed;
   }
   Honeycomb::logger(Honeycomb::Logger::INFO,
                     std::format("  Total time elapsed: {:.4e} (s)", total_time * 1.0e-3));

   std::FILE *fp_2 = std::fopen(file_name.c_str(), "w");

   for (size_t i = 0; i < xBjs.size(); i++) {
      double res = 0;
      for (size_t j = 1; j < sol_fin._distr_p.size(); j++) {
         double pref = 1.0 / 9.0;
         if (j % 2 == 0) pref *= 4; // Up-type quarks are j=2, 4, 6
         res += pref * sol_fin._distr_p[j].dot(weights[i].weights);
      }

      std::fprintf(fp_2, "%.16e\t%.16e\n", xBjs[i], res);
   }
   std::fclose(fp_2);
   Honeycomb::logger(Honeycomb::Logger::INFO,
                     std::format("================== End: {:s} =================", file_name));
}