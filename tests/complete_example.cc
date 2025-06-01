#include <honeycomb2/honeycomb2.hpp>
#include "alpha_s.hpp"
#include "honeycomb_120_25_grid_points.hpp"
#include "original_model_functions.hpp"
#include "solution.hpp"
#include "utilities.hpp"

int main()
{
   // For timing purposes
   Honeycomb::timer::mark begin = Honeycomb::timer::now();
   Honeycomb::timer::mark end   = Honeycomb::timer::now();

   // Grid setup
   const size_t n    = 7;
   const double rmin = 0.001;
   Honeycomb::Grid2D grid;
   grid = Honeycomb::generate_compliant_Grid2D(n, {rmin, 0.1, 0.4, 1}, {12, 8, 7});
   Honeycomb::Discretization discr(grid);

   // Print grid information
   Honeycomb::logger(Honeycomb::Logger::INFO,
                     std::format("Total grid size: {:d}x{:d}", grid.size, grid.size));
   Honeycomb::logger(Honeycomb::Logger::INFO,
                     std::format("Total x123 size: {:d}x{:d}", grid.c_size, grid.c_size));

   // Load/compute evolution kernels, timing it
   Honeycomb::logger(Honeycomb::Logger::INFO, std::format("Kernel Discretization..."));
   begin                   = Honeycomb::timer::now();
   const double Nc         = 3; // NC = 1 for tests
   Honeycomb::Kernels kers = Honeycomb::load_kernels("full_kernels.cereal", grid, Nc);
   end                     = Honeycomb::timer::now();
   Honeycomb::logger(Honeycomb::Logger::INFO,
                     std::format("  Elapsed: {:.4e} (ms)", Honeycomb::timer::elapsed_ms(end, begin)));

   // Initial model setup
   Honeycomb::InputModel model;

   model.SetModel(Honeycomb::InputModel::T_UP, orignal_models::Tu_test);
   model.SetModel(Honeycomb::InputModel::DT_UP, orignal_models::DTu_test);

   model.SetModel(Honeycomb::InputModel::T_DN, orignal_models::Td_test);
   model.SetModel(Honeycomb::InputModel::DT_DN, orignal_models::DTd_test);

   model.SetModel(Honeycomb::InputModel::T_ST, orignal_models::Ts_test);
   model.SetModel(Honeycomb::InputModel::DT_ST, orignal_models::DTs_test);

   model.SetModel(Honeycomb::InputModel::T_P_GL, orignal_models::TFp_test);
   model.SetModel(Honeycomb::InputModel::T_M_GL, orignal_models::TFm_test);

   // Evolution prefactor, DO NOT CHANGE
   const double pref = -1.0;

   // Initial scale, GeV^2
   const double Q02 = 1.0;
   const double t0  = log(Q02);

   // Final scale, GeV^2
   const double Qf2 = 1.0e+4;

   // Mass^2 of quarks (down, up, strange, charm, bottom, top)
   std::array<double, 6> thresholds = {0, 0, 0, 1.6129, 17.4724, 1.0e+6};

   // Setup for evolution
   // \alpha_s / 4\pi
   auto as = Honeycomb::GetAlphaS_o_4pi(thresholds);

   // Generate 'Solution' at initial scale and list of intermediate scales for evolution
   auto [inter_scale, sol1] = get_initial_solution(Q02, Qf2, thresholds, &discr, model);

   // Copy only the log(m_q^2) scales, to push new flavors in evolution
   std::vector<double> log_active_thresholds(inter_scale.begin(), inter_scale.end() - 1);

   std::vector<std::pair<double, Honeycomb::Solution>> solutions;

   // Setup callback for evolution to add new flavor at the threshold
   auto callback = [&log_active_thresholds, &solutions](double t, const Honeycomb::Kernels &,
                                                        Honeycomb::Solution &S) -> void {
      auto it = std::find_if(log_active_thresholds.begin(), log_active_thresholds.end(), [t](double x) {
         return std::abs(x - t) < 1.0e-14;
      });

      // t = log(m_q^2)
      if (it != log_active_thresholds.end()) S.PushFlavor();

      solutions.emplace_back(std::make_pair(t, S));
      return;
   };

   // Construct the evolver, object which will perform the actual evolution
   Honeycomb::runge_kutta::GenericRungeKutta<Honeycomb::Kernels, Honeycomb::Solution, 13> evolver(
       kers, sol1, Honeycomb::runge_kutta::DOPRI8, as, pref, t0, 0.01, callback);

   // Evolve from initial to final scale, between each threshold uses 40 steps
   Honeycomb::logger(Honeycomb::Logger::INFO, std::format("Evolving..."));
   begin = Honeycomb::timer::now();
   evolver(inter_scale, 40);
   end = Honeycomb::timer::now();
   Honeycomb::logger(Honeycomb::Logger::INFO,
                     std::format("  Elapsed: {:.4e} (ms)", Honeycomb::timer::elapsed_ms(end, begin)));

   // Extract solution at final scale from the Evolver
   Honeycomb::Solution sol_fin = evolver.GetSolution();

   // The solution extracted from the evolver must be equal to the last
   // entry of the solutions array to which we appendend in the callback function.
   // This is a sanity check.
   bool are_equal = sol_fin.is_equalt_to(solutions.back().second);
   if (are_equal) {
      Honeycomb::logger(Honeycomb::Logger::INFO, "Solutions match!");
   } else {

      Honeycomb::logger(Honeycomb::Logger::WARNING, "Solutions do not match!");
   }

   // Put solution into physical basis
   sol1.RotateToPhysicalBasis();
   sol_fin.RotateToPhysicalBasis();

   // Obtain weights for 'd2' observable
   // This needs to be convoluted with S^+, not T!
   Honeycomb::logger(Honeycomb::Logger::INFO, std::format("d2 Weights computation"));
   begin = Honeycomb::timer::now();
   Honeycomb::D2Weights d2_weights(grid);
   end = Honeycomb::timer::now();
   Honeycomb::logger(Honeycomb::Logger::INFO,
                     std::format("  Elapsed: {:.4e} (ms)", Honeycomb::timer::elapsed_ms(end, begin)));

   // Compute d2 for each flavor and print it to file
   {
      // Open file in "w" = write mode: overwrites old content if any
      std::FILE *fp = std::fopen("d2_at_1e+4_GeV2.dat", "w");
      // Single column, rows are, in order: Down, Up, Strange, ... until last non-zero flavor
      for (size_t i = 1; i < sol_fin._distr_p.size(); i++) {
         std::fprintf(fp, "%.16e\n", d2_weights.ComputeSingleQuark(sol_fin._distr_p[i]));
      }
      std::fclose(fp);
   }

   // Obtain weights for 'd1' observable -> Efremov-LeaderTeryaev sum-rule (up to factor 2)
   // This needs to be convoluted with S^+, not T!
   Honeycomb::logger(Honeycomb::Logger::INFO, std::format("d1 Weights computation"));
   begin = Honeycomb::timer::now();
   Honeycomb::ELTWeights elt_weights(grid);
   end = Honeycomb::timer::now();
   Honeycomb::logger(Honeycomb::Logger::INFO,
                     std::format("  Elapsed: {:.4e} (ms)", Honeycomb::timer::elapsed_ms(end, begin)));

   // Compute d1 for each flavor and print it to file
   {
      // Open file in "w" = write mode: overwrites old content if any
      std::FILE *fp = std::fopen("d1_at_1e+4_GeV2.dat", "w");
      // Single column, rows are, in order: Down, Up, Strange, ... until last non-zero flavor
      for (size_t i = 1; i < sol_fin._distr_p.size(); i++) {
         std::fprintf(fp, "%.16e\n", elt_weights.ComputeSingleQuark(sol_fin._distr_p[i]));
      }
      std::fclose(fp);
   }

   // Save the initial and evolved solutions to file in the T, DT basis
   {
      // Construct T and DT from S^{\pm}
      Honeycomb::OutputModel out_model(sol_fin);
      Honeycomb::OutputModel out_model_ini(sol1);

      std::FILE *fp   = std::fopen("T_and_DT_final_scale.dat", "w");
      std::FILE *fp_1 = std::fopen("T_and_DT_initial_scale.dat", "w");

      // Get the points from 'original' honeycomb grid
      for (Honeycomb::RnC::Triplet &x123 : honeycomb_points) {

         // Column 1, 2, and 3 are x1, x2, x3
         std::fprintf(fp, "%.16e\t%.16e\t%.16e\t", x123(0), x123(1), x123(2));

         // Save distributions. Order is:
         // T_{3F}^+, T_{3F}^-, T_d, DT_d, T_u, DT_u, ...
         // ^^^^^^^^^^^^^^^^^^
         // These are in Eq. 23 of 2405.01162

         //------------------------------------------
         // This is final scale model
         std::fprintf(fp, "%.16e\t", out_model.GetDistribution(Honeycomb::OutputModel::T_P_GL, x123));
         std::fprintf(fp, "%.16e\t", out_model.GetDistribution(Honeycomb::OutputModel::T_M_GL, x123));

         std::fprintf(fp, "%.16e\t", out_model.GetDistribution(Honeycomb::OutputModel::T_DN, x123));
         std::fprintf(fp, "%.16e\t", out_model.GetDistribution(Honeycomb::OutputModel::T_DN, x123));

         std::fprintf(fp, "%.16e\t", out_model.GetDistribution(Honeycomb::OutputModel::T_UP, x123));
         std::fprintf(fp, "%.16e\t", out_model.GetDistribution(Honeycomb::OutputModel::T_UP, x123));

         std::fprintf(fp, "%.16e\t", out_model.GetDistribution(Honeycomb::OutputModel::T_ST, x123));
         std::fprintf(fp, "%.16e\t", out_model.GetDistribution(Honeycomb::OutputModel::T_ST, x123));

         std::fprintf(fp, "%.16e\t", out_model.GetDistribution(Honeycomb::OutputModel::T_CH, x123));
         std::fprintf(fp, "%.16e\t", out_model.GetDistribution(Honeycomb::OutputModel::T_CH, x123));

         std::fprintf(fp, "%.16e\t", out_model.GetDistribution(Honeycomb::OutputModel::T_BM, x123));
         std::fprintf(fp, "%.16e\t", out_model.GetDistribution(Honeycomb::OutputModel::T_BM, x123));

         // Omit Top quark, essentially irrelevant for us.
         std::fprintf(fp, "\n");

         //------------------------------------------
         // This is initial scale model

         std::fprintf(fp_1, "%.16e\t", out_model_ini.GetDistribution(Honeycomb::OutputModel::T_P_GL, x123));
         std::fprintf(fp_1, "%.16e\t", out_model_ini.GetDistribution(Honeycomb::OutputModel::T_M_GL, x123));
         std::fprintf(fp_1, "%.16e\t", out_model_ini.GetDistribution(Honeycomb::OutputModel::T_DN, x123));
         std::fprintf(fp_1, "%.16e\t", out_model_ini.GetDistribution(Honeycomb::OutputModel::T_DN, x123));
         std::fprintf(fp_1, "%.16e\t", out_model_ini.GetDistribution(Honeycomb::OutputModel::T_UP, x123));
         std::fprintf(fp_1, "%.16e\t", out_model_ini.GetDistribution(Honeycomb::OutputModel::T_UP, x123));
         std::fprintf(fp_1, "%.16e\t", out_model_ini.GetDistribution(Honeycomb::OutputModel::T_ST, x123));
         std::fprintf(fp_1, "%.16e\t", out_model_ini.GetDistribution(Honeycomb::OutputModel::T_ST, x123));
         std::fprintf(fp_1, "%.16e\t", out_model_ini.GetDistribution(Honeycomb::OutputModel::T_CH, x123));
         std::fprintf(fp_1, "%.16e\t", out_model_ini.GetDistribution(Honeycomb::OutputModel::T_CH, x123));
         std::fprintf(fp_1, "%.16e\t", out_model_ini.GetDistribution(Honeycomb::OutputModel::T_BM, x123));
         std::fprintf(fp_1, "%.16e\t", out_model_ini.GetDistribution(Honeycomb::OutputModel::T_BM, x123));

         // Omit Top quark, essentially irrelevant for us.
         std::fprintf(fp_1, "\n");
      }
      std::fclose(fp);
      std::fclose(fp_1);
   }

   return 0;
}