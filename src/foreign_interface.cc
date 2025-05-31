#include "discretization.hpp"
#include "obs_weights.hpp"
#include "timer.hpp"
#include "utilities.hpp"
#include <honeycomb2/foreign_interface.hpp>
#include <thread>

static Honeycomb::ForeignInterfaceState state;

namespace Honeycomb
{

bool _compare_thr(const std::vector<double> &v1, const std::vector<double> &v2)
{
   bool res = true;
   for (size_t i = 0; i < 6; i++) {
      res = res && is_near(v1[i], v2[i], 1.0e-12);
   }
   return res;
}

void set_up(const std::string &config_name)
{
   Honeycomb::timer::mark begin = Honeycomb::timer::now();
   Honeycomb::timer::mark end   = Honeycomb::timer::now();

   logger(Logger::INFO, "Parsing config file...", false);
   begin = timer::now();
   ConfigParser cp(read_file_to_str(config_name));
   end = timer::now();
   logger(Logger::INFO, std::format("  Elapsed: {:.4e} (ms)", timer::elapsed_ms(end, begin)));

   // Parse the grid info
   const size_t n               = cp.GetValue<size_t>("n_phi");
   const double rmin            = cp.GetValue<double>("rmin");
   std::vector<double> r_values = cp.GetValue<std::vector<double>>("r_v");
   std::vector<size_t> r_interv = cp.GetValue<std::vector<size_t>>("r_i");
   r_values.insert(r_values.begin(), rmin);

   state.grid = generate_compliant_Grid2D(n, r_values, r_interv);
   logger(Logger::INFO,
          std::format("Total grid size: {:d}x{:d}", state.discr._grid.size, state.discr._grid.size));
   logger(Logger::INFO,
          std::format("Total x123 size: {:d}x{:d}", state.discr._grid.c_size, state.discr._grid.c_size));

   if (!load_weights(cp.GetValue("wf", 0), state.d2_weights)) {
      state.d2_weights.GetWeights();
      save_weights(state.d2_weights, cp.GetValue("wf", 0));
   }
   if (!load_weights(cp.GetValue("wf", 1), state.elt_weights)) {
      state.elt_weights.GetWeights();
      save_weights(state.elt_weights, cp.GetValue("wf", 1));
   }

   state.interm_scales = cp.GetValue<std::vector<double>>("Scales");
   state.thresholds    = cp.GetValue<std::vector<double>>("mthr");
   if (state.thresholds.size() != 6) {
      logger(Logger::ERROR, std::format("Config must contain exactly 6 mass thresholds. {:d} found.",
                                        state.thresholds.size()));
   }

   std::vector<std::string> eo_file_list = cp.GetMap().at("eo");
   if (eo_file_list.size() != state.interm_scales.size() - 1) {
      logger(Logger::ERROR, std::format("Mismatch between number of finale scales ({:d}) requested and list "
                                        "of files for evolution operators ({:d}).",
                                        state.interm_scales.size() - 1, eo_file_list.size()));
   }

   std::string ker_file = "fi_kernels.cereal";

   if (cp.GetMap().find("kf") != cp.GetMap().end()) {
      ker_file = cp.GetMap().at("kf")[0];
   }

   for (size_t i = 0; i < eo_file_list.size(); i++) {
      std::string f_name          = eo_file_list[i];
      auto [ok_load, eo_load_tmp] = load_evolution_operator(f_name, &state.grid);

      const double Q02 = state.interm_scales[i];
      const double Qf2 = state.interm_scales[i + 1];
      const double t0  = log(Q02);
      const double tF  = log(Qf2);

      // Did not load properly OR scales do not match
      if (!ok_load || !is_near(t0, eo_load_tmp.t0tF.first, 1.0e-10)
          || !is_near(tF, eo_load_tmp.t0tF.second, 1.0e-10)
          || !_compare_thr(eo_load_tmp._thresholds, state.thresholds)) {
         // logger(Logger::ERROR, "TODO");
         Kernels kers = Honeycomb::load_kernels(ker_file, state.grid, 3);
         std::array<double, 6> tmp_thr;
         for (size_t j = 0; j < 6; j++) {
            tmp_thr[j] = state.thresholds[j];
         }

         logger(Logger::INFO, "Computing evolution operator...");
         begin       = timer::now();
         eo_load_tmp = compute_evolution_operator(&state.grid, kers, Q02, Qf2, tmp_thr, state.as_fnc);
         save_evolution_operator(eo_load_tmp, f_name);
         end = timer::now();
         logger(Logger::INFO, std::format("  Elapsed: {:.4e} (ms)", timer::elapsed_ms(end, begin)));
      }
      state.evol_op.emplace_back(eo_load_tmp);
   }
}
}; // namespace Honeycomb

extern "C" void set_up_(const char *config_name, int len)
{
   (void)len;
   Honeycomb::timer::mark begin = Honeycomb::timer::now();
   Honeycomb::set_up(std::string(config_name));
   Honeycomb::timer::mark end = Honeycomb::timer::now();
   Honeycomb::logger(Honeycomb::Logger::INFO,
                     std::format("Full setup took {:.4e} (ms)", Honeycomb::timer::elapsed_ms(end, begin)));
   std::this_thread::sleep_for(std::chrono::seconds(30));
}