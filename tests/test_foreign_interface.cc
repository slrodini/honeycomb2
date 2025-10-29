#include <honeycomb2/honeycomb2.hpp>
#include <honeycomb2/honeycomb2_c_api.h>
#include <test_config.hpp>

int check_initial_condition();

int main()
{
   std::string file = TESTS_PATH "/example.config";

   hc2_fi_set_up_(file.c_str(), file.size());
   hc2_fi_evolve_();
   hc2_fi_unload_();

   // Grid setup
   const size_t n    = 7;
   const double rmin = 0.001;
   Honeycomb::Grid2D grid;
   grid = Honeycomb::generate_compliant_Grid2D(n, {rmin, 0.1, 0.4, 1}, {12, 8, 7});
   Honeycomb::Discretization discr(grid);

   // Load/compute evolution kernels, timing it
   const double Nc         = 3; // NC = 1 for tests
   Honeycomb::Kernels kers = Honeycomb::load_kernels("fi_kernels.cereal", grid, Nc);

   // Initial model setup
   Honeycomb::InputModel model = Honeycomb::PreImplementedModels::GetModel("pim_LFWA_fitted");

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
   std::vector<double> log_active_thresholds(inter_scale.begin() + 1, inter_scale.end() - 1);

   inter_scale.push_back(log(4));
   inter_scale.push_back(log(25));
   inter_scale.push_back(log(100));
   std::sort(inter_scale.begin(), inter_scale.end());

   auto tableau = Honeycomb::rk::PreImplementedTableau::DOPRI8;

   auto rhs = Honeycomb::Get_Solution_rk_rhs(as, kers);

   using TInfo = Honeycomb::rk::TimeInfo;
   // Construct the evolver, object which will perform the actual evolution
   Honeycomb::rk::RungeKutta<Honeycomb::Solution, tableau.stages> evolver(
       tableau, sol1, TInfo{inter_scale, 0.1, 10}, rhs);

   std::vector<std::pair<double, Honeycomb::Solution>> solutions;

   // Setup callback for evolution to add new flavor at the threshold
   std::function<void(double, const TInfo &, Honeycomb::Solution &)> callback
       = [&log_active_thresholds, &solutions](double t, const TInfo &,
                                              Honeycomb::Solution &S) -> void {
      auto it
          = std::find_if(log_active_thresholds.begin(), log_active_thresholds.end(), [t](double x) {
               return std::abs(x - t) < 1.0e-14;
            });

      // t = log(m_q^2)
      if (it != log_active_thresholds.end()) S.PushFlavor();
      else solutions.emplace_back(std::make_pair(t, S));
      return;
   };

   evolver.AddCallback(callback);

   // Evolve from initial to final scale, between each threshold uses 40 steps
   evolver();

   // Extract solution at final scale from the Evolver
   Honeycomb::Solution sol_fin = evolver.GetSolution();

   // The solution extracted from the evolver must be equal to the last
   // entry of the solutions array to which we appendend in the callback function.
   // This is a sanity check.

   Honeycomb::logger(Honeycomb::Logger::INFO,
                     "Number of cached solution is: " + std::to_string(solutions.size()));

   for (size_t i = 0; i < solutions.size(); i++) {
      solutions[i].second.RotateToPhysicalBasis();
   }

   std::vector<double> check_scales = {4, 25, 100, 10000};
   std::vector<Honeycomb::OutputModel> out_models;
   for (size_t i = 0; i < solutions.size(); i++) {
      if (!Honeycomb::is_near(log(check_scales[i]), solutions[i].first, 1.0e-10)) {
         Honeycomb::logger(Honeycomb::Logger::ERROR,
                           std::format("Mismatch in the scales. {:.10e}, {:.10e}",
                                       log(check_scales[i]), solutions[i].first));
      }
      out_models.emplace_back(Honeycomb::OutputModel(solutions[i].second));
   }

   double m = 0;
   for (size_t k = 0; k < check_scales.size(); k++) {
      for (int i = 0; i < 13; i++) {
         for (size_t j = 0; j < grid._x123.size(); j++) {
            double x1 = grid._x123[j][0];
            double x2 = grid._x123[j][1];
            double x3 = grid._x123[j][2];
            double fi = hc2_fi_get_model_(&i, &check_scales[k], &x1, &x2, &x3);

            Honeycomb::OutputModel::FNC d = static_cast<Honeycomb::OutputModel::FNC>(i);

            double native = out_models[k].GetDistribution(d, x1, x2, x3);

            m = std::max(m, std::fabs(fi - native));
         }
      }
   }
   if (!Honeycomb::is_near(m, 0., 1.0e-12)) {
      Honeycomb::logger(Honeycomb::Logger::ERROR,
                        std::format("Max difference is too large: {:.10e}", m));
   }

   return 0;
}
