#include <honeycomb2/honeycomb2.hpp>

#include "honeycomb_120_25_grid_points.hpp"

#include "solution.hpp"
#include "test_rd_T_and_DT_final_scale.hpp"
#include "test_rd_T_and_DT_initial_scale.hpp"

static int ret_val = 0;

#define GEN_REGRESSION 0

void generate_regression(const Honeycomb::OutputModel &final,
                         const Honeycomb::OutputModel &initial);

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

   // Load/compute evolution kernels, timing it
   Honeycomb::logger(Honeycomb::Logger::INFO, std::format("Kernel Discretization..."));
   begin                   = Honeycomb::timer::now();
   const double Nc         = 3; // NC = 1 for tests
   Honeycomb::Kernels kers = Honeycomb::load_kernels("full_kernels.cereal", grid, Nc);
   end                     = Honeycomb::timer::now();
   Honeycomb::logger(
       Honeycomb::Logger::INFO,
       std::format("  Elapsed: {:.4e} (ms)", Honeycomb::timer::elapsed_ms(end, begin)));

   // Initial model setup
   Honeycomb::InputModel model = Honeycomb::PreImplementedModels::GetModel("pim_original");

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
      auto it
          = std::find_if(log_active_thresholds.begin(), log_active_thresholds.end(), [t](double x) {
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
   Honeycomb::logger(
       Honeycomb::Logger::INFO,
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
      ret_val += 1;
   }
   if (ret_val != 0) return ret_val;

   // Put solution into physical basis
   sol1.RotateToPhysicalBasis();
   sol_fin.RotateToPhysicalBasis();

   // Save the initial and evolved solutions to file in the T, DT basis
   {

      Honeycomb::OutputModel final(sol_fin);
      Honeycomb::OutputModel initial(sol1);
      if (GEN_REGRESSION == 1) {
         generate_regression(final, initial);
         return 0;
      } else {
         size_t index = 0;
         if (honeycomb_points.size() != T_DT_regression_final.size()
             || honeycomb_points.size() != T_DT_regression_initial.size()) {
            Honeycomb::logger(Honeycomb::Logger::ERROR, "Regression data have wrong shape.");
         }
         for (Honeycomb::RnC::Triplet &x123 : honeycomb_points) {
            double curr_pt = 0;
            curr_pt        = final.GetDistribution(Honeycomb::OutputModel::T_P_GL, x123);
            if (!Honeycomb::is_near(curr_pt, T_DT_regression_final[index][3], 1.0e-6)) ret_val += 1;
            curr_pt = final.GetDistribution(Honeycomb::OutputModel::T_M_GL, x123);
            if (!Honeycomb::is_near(curr_pt, T_DT_regression_final[index][4], 1.0e-6)) ret_val += 1;
            curr_pt = final.GetDistribution(Honeycomb::OutputModel::T_DN, x123);
            if (!Honeycomb::is_near(curr_pt, T_DT_regression_final[index][5], 1.0e-6)) ret_val += 1;
            curr_pt = final.GetDistribution(Honeycomb::OutputModel::DT_DN, x123);
            if (!Honeycomb::is_near(curr_pt, T_DT_regression_final[index][6], 1.0e-6)) ret_val += 1;
            curr_pt = final.GetDistribution(Honeycomb::OutputModel::T_UP, x123);
            if (!Honeycomb::is_near(curr_pt, T_DT_regression_final[index][7], 1.0e-6)) ret_val += 1;
            curr_pt = final.GetDistribution(Honeycomb::OutputModel::DT_UP, x123);
            if (!Honeycomb::is_near(curr_pt, T_DT_regression_final[index][8], 1.0e-6)) ret_val += 1;
            curr_pt = final.GetDistribution(Honeycomb::OutputModel::T_ST, x123);
            if (!Honeycomb::is_near(curr_pt, T_DT_regression_final[index][9], 1.0e-6)) ret_val += 1;
            curr_pt = final.GetDistribution(Honeycomb::OutputModel::DT_ST, x123);
            if (!Honeycomb::is_near(curr_pt, T_DT_regression_final[index][10], 1.0e-6))
               ret_val += 1;
            curr_pt = final.GetDistribution(Honeycomb::OutputModel::T_CH, x123);
            if (!Honeycomb::is_near(curr_pt, T_DT_regression_final[index][11], 1.0e-6))
               ret_val += 1;
            curr_pt = final.GetDistribution(Honeycomb::OutputModel::DT_CH, x123);
            if (!Honeycomb::is_near(curr_pt, T_DT_regression_final[index][12], 1.0e-6))
               ret_val += 1;
            curr_pt = final.GetDistribution(Honeycomb::OutputModel::T_BM, x123);
            if (!Honeycomb::is_near(curr_pt, T_DT_regression_final[index][13], 1.0e-6))
               ret_val += 1;
            curr_pt = final.GetDistribution(Honeycomb::OutputModel::DT_BM, x123);
            if (!Honeycomb::is_near(curr_pt, T_DT_regression_final[index][14], 1.0e-6))
               ret_val += 1;
            curr_pt = initial.GetDistribution(Honeycomb::OutputModel::T_P_GL, x123);
            if (!Honeycomb::is_near(curr_pt, T_DT_regression_initial[index][3], 1.0e-6))
               ret_val += 1;
            curr_pt = initial.GetDistribution(Honeycomb::OutputModel::T_M_GL, x123);
            if (!Honeycomb::is_near(curr_pt, T_DT_regression_initial[index][4], 1.0e-6))
               ret_val += 1;
            curr_pt = initial.GetDistribution(Honeycomb::OutputModel::T_DN, x123);
            if (!Honeycomb::is_near(curr_pt, T_DT_regression_initial[index][5], 1.0e-6))
               ret_val += 1;
            curr_pt = initial.GetDistribution(Honeycomb::OutputModel::DT_DN, x123);
            if (!Honeycomb::is_near(curr_pt, T_DT_regression_initial[index][6], 1.0e-6))
               ret_val += 1;
            curr_pt = initial.GetDistribution(Honeycomb::OutputModel::T_UP, x123);
            if (!Honeycomb::is_near(curr_pt, T_DT_regression_initial[index][7], 1.0e-6))
               ret_val += 1;
            curr_pt = initial.GetDistribution(Honeycomb::OutputModel::DT_UP, x123);
            if (!Honeycomb::is_near(curr_pt, T_DT_regression_initial[index][8], 1.0e-6))
               ret_val += 1;
            curr_pt = initial.GetDistribution(Honeycomb::OutputModel::T_ST, x123);
            if (!Honeycomb::is_near(curr_pt, T_DT_regression_initial[index][9], 1.0e-6))
               ret_val += 1;
            curr_pt = initial.GetDistribution(Honeycomb::OutputModel::DT_ST, x123);
            if (!Honeycomb::is_near(curr_pt, T_DT_regression_initial[index][10], 1.0e-6))
               ret_val += 1;
            curr_pt = initial.GetDistribution(Honeycomb::OutputModel::T_CH, x123);
            if (!Honeycomb::is_near(curr_pt, T_DT_regression_initial[index][11], 1.0e-6))
               ret_val += 1;
            curr_pt = initial.GetDistribution(Honeycomb::OutputModel::DT_CH, x123);
            if (!Honeycomb::is_near(curr_pt, T_DT_regression_initial[index][12], 1.0e-6))
               ret_val += 1;
            curr_pt = initial.GetDistribution(Honeycomb::OutputModel::T_BM, x123);
            if (!Honeycomb::is_near(curr_pt, T_DT_regression_initial[index][13], 1.0e-6))
               ret_val += 1;
            curr_pt = initial.GetDistribution(Honeycomb::OutputModel::DT_BM, x123);
            if (!Honeycomb::is_near(curr_pt, T_DT_regression_initial[index][14], 1.0e-6))
               ret_val += 1;

            if (ret_val != 0) return ret_val;
            index++;
         }
      }
   }

   return 0;
}

void generate_regression(const Honeycomb::OutputModel &final, const Honeycomb::OutputModel &initial)
{
   std::FILE *fp   = std::fopen("T_and_DT_final_scale.hpp", "w");
   std::FILE *fp_1 = std::fopen("T_and_DT_initial_scale.hpp", "w");

   std::fprintf(fp, "#include <vector>\n\n"
                    "std::vector<std::vector<double>> T_DT_regression_final = {\n");
   std::fprintf(fp_1, "#include <vector>\n\n"
                      "std::vector<std::vector<double>> T_DT_regression_initial = {\n");

   // Get the points from 'original' honeycomb grid
   for (Honeycomb::RnC::Triplet &x123 : honeycomb_points) {

      // Column 1, 2, and 3 are x1, x2, x3
      std::fprintf(fp, "   {%.16e, %.16e, %.16e, ", x123(0), x123(1), x123(2));
      std::fprintf(fp_1, "   {%.16e, %.16e, %.16e, ", x123(0), x123(1), x123(2));

      // Save distributions. Order is:
      // T_{3F}^+, T_{3F}^-, T_d, DT_d, T_u, DT_u, ...
      // ^^^^^^^^^^^^^^^^^^
      // These are in Eq. 23 of 2405.01162

      //------------------------------------------
      // This is final scale model
      std::fprintf(fp, "%.16e, ", final.GetDistribution(Honeycomb::OutputModel::T_P_GL, x123));
      std::fprintf(fp, "%.16e, ", final.GetDistribution(Honeycomb::OutputModel::T_M_GL, x123));

      std::fprintf(fp, "%.16e, ", final.GetDistribution(Honeycomb::OutputModel::T_DN, x123));
      std::fprintf(fp, "%.16e, ", final.GetDistribution(Honeycomb::OutputModel::DT_DN, x123));

      std::fprintf(fp, "%.16e, ", final.GetDistribution(Honeycomb::OutputModel::T_UP, x123));
      std::fprintf(fp, "%.16e, ", final.GetDistribution(Honeycomb::OutputModel::DT_UP, x123));

      std::fprintf(fp, "%.16e, ", final.GetDistribution(Honeycomb::OutputModel::T_ST, x123));
      std::fprintf(fp, "%.16e, ", final.GetDistribution(Honeycomb::OutputModel::DT_ST, x123));

      std::fprintf(fp, "%.16e, ", final.GetDistribution(Honeycomb::OutputModel::T_CH, x123));
      std::fprintf(fp, "%.16e, ", final.GetDistribution(Honeycomb::OutputModel::DT_CH, x123));

      std::fprintf(fp, "%.16e, ", final.GetDistribution(Honeycomb::OutputModel::T_BM, x123));
      std::fprintf(fp, "%.16e, ", final.GetDistribution(Honeycomb::OutputModel::DT_BM, x123));

      // Omit Top quark, essentially irrelevant for us.
      std::fprintf(fp, "},\n");

      //------------------------------------------
      // This is initial scale model

      std::fprintf(fp_1, "%.16e, ", initial.GetDistribution(Honeycomb::OutputModel::T_P_GL, x123));
      std::fprintf(fp_1, "%.16e, ", initial.GetDistribution(Honeycomb::OutputModel::T_M_GL, x123));
      std::fprintf(fp_1, "%.16e, ", initial.GetDistribution(Honeycomb::OutputModel::T_DN, x123));
      std::fprintf(fp_1, "%.16e, ", initial.GetDistribution(Honeycomb::OutputModel::DT_DN, x123));
      std::fprintf(fp_1, "%.16e, ", initial.GetDistribution(Honeycomb::OutputModel::T_UP, x123));
      std::fprintf(fp_1, "%.16e, ", initial.GetDistribution(Honeycomb::OutputModel::DT_UP, x123));
      std::fprintf(fp_1, "%.16e, ", initial.GetDistribution(Honeycomb::OutputModel::T_ST, x123));
      std::fprintf(fp_1, "%.16e, ", initial.GetDistribution(Honeycomb::OutputModel::DT_ST, x123));
      std::fprintf(fp_1, "%.16e, ", initial.GetDistribution(Honeycomb::OutputModel::T_CH, x123));
      std::fprintf(fp_1, "%.16e, ", initial.GetDistribution(Honeycomb::OutputModel::DT_CH, x123));
      std::fprintf(fp_1, "%.16e, ", initial.GetDistribution(Honeycomb::OutputModel::T_BM, x123));
      std::fprintf(fp_1, "%.16e, ", initial.GetDistribution(Honeycomb::OutputModel::DT_BM, x123));

      // Omit Top quark, essentially irrelevant for us.
      std::fprintf(fp_1, "},\n");
   }
   std::fprintf(fp, "};\n");
   std::fprintf(fp_1, "};\n");

   std::fclose(fp);
   std::fclose(fp_1);
}