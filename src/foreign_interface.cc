#include <honeycomb2/foreign_interface.hpp>
#include <honeycomb2/alpha_s.hpp>
#include <honeycomb2/honeycomb2_c_api.h>

static Honeycomb::ForeignInterfaceState state;

namespace Honeycomb
{

void set_model(const std::string &what)
{
   bool do_I_have_it = Honeycomb::PreImplementedModels::QueryModelAvailability(what);
   if (do_I_have_it) {
      state._models = Honeycomb::PreImplementedModels::GetModel(what);
      Honeycomb::logger(Honeycomb::Logger::INFO, "ForeignInterface: using model <" + what + ">.");
   } else {
      Honeycomb::logger(Honeycomb::Logger::ERROR,
                        "ForeignInterface: model <" + what
                            + "> is not available in the list of pre-implemented models.");
   }
}

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
   timer::mark begin = timer::now();
   timer::mark end   = timer::now();

   logger(Logger::INFO, "Parsing config file...");
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

   state.discr = new Discretization(state.grid);

   logger(Logger::INFO, std::format("Total grid size: {:d}x{:d}", state.discr->_grid.size,
                                    state.discr->_grid.size));
   logger(Logger::INFO, std::format("Total x123 size: {:d}x{:d}", state.discr->_grid.c_size,
                                    state.discr->_grid.c_size));

   const auto entries = cp.GetMap();

   if (entries.find("d2f") != entries.end()) {
      state.d2_weights = new D2Weights(state.grid, false);
      if (!load_weights(cp.GetValue("d2f", 0), *state.d2_weights)) {
         state.d2_weights->GetWeights();
         save_weights(*state.d2_weights, cp.GetValue("d2f", 0));
      }
   }

   if (entries.find("eltf") != entries.end()) {
      state.elt_weights = new ELTWeights(state.grid, false);

      if (!load_weights(cp.GetValue("eltf", 0), *state.elt_weights)) {
         state.elt_weights->GetWeights();
         save_weights(*state.elt_weights, cp.GetValue("eltf", 0));
      }
   }

   if (entries.find("g2f") != entries.end()) {
      state.v_g2_weights = new VectorG2Weights(state.grid, false);

      if (!load_weights(cp.GetValue("g2f", 0), *state.v_g2_weights)) {
         state.v_g2_weights->GetWeights();
         save_weights(*state.v_g2_weights, cp.GetValue("g2f", 0));
      }
   }

   state.interm_scales = cp.GetValue<std::vector<double>>("Scales");
   if (state.interm_scales.size() <= 1) {
      logger(Logger::ERROR, "Not enough scales to evolve. At least two scales "
                            "(initial and final) must be specified.");
   }

   state.thresholds = cp.GetValue<std::vector<double>>("mthr");
   if (state.thresholds.size() != 6) {
      logger(Logger::ERROR,
             std::format("Config must contain exactly 6 mass thresholds. {:d} found.",
                         state.thresholds.size()));
   }

   // Compute initial nf
   for (size_t i = 0; i < 5; i++) {
      if (state.thresholds[i + 1] < state.thresholds[i]) {
         logger(Logger::ERROR,
                "Thresholds are not increasing, backward evolution not supported yet.");
      }
   }

   std::array<double, 6> tmp_thr;
   for (size_t j = 0; j < 6; j++) {
      tmp_thr[j] = state.thresholds[j];
   }
   state.as_fnc = GetAlphaS_o_4pi(tmp_thr);

   state.nf_initial_scale = 0;
   for (size_t i = 0; i < state.thresholds.size(); i++) {
      if (state.thresholds[i] <= state.interm_scales[0]) {
         state.nf_initial_scale++;
      } else {
         break;
      }
   }

   std::vector<std::string> eo_file_list = cp.GetMap().at("eo");
   if (eo_file_list.size() != state.interm_scales.size() - 1) {
      logger(Logger::ERROR,
             std::format("Mismatch between number of finale scales ({:d}) requested and list "
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

         logger(Logger::INFO, "Computing evolution operator...");
         begin = timer::now();
         eo_load_tmp
             = compute_evolution_operator(&state.grid, kers, Q02, Qf2, tmp_thr, state.as_fnc);
         save_evolution_operator(eo_load_tmp, f_name);
         end = timer::now();
         logger(Logger::INFO, std::format("  Elapsed: {:.4e} (ms)", timer::elapsed_ms(end, begin)));
      }
      state.evol_op.emplace_back(eo_load_tmp);
   }

   if (cp.GetMap().find("model") != cp.GetMap().end()) {
      std::vector<std::string> models = cp.GetMap().at("model");
      if (models.size() > 1) {
         logger(Logger::INFO, "Mulitple models specified, but interface support only evolution of "
                              "single model. Using the last entry.");
      }
      std::string model = models.back();
      set_model(model);

   } else {
      logger(Logger::INFO,
             "No pre-implemented model has been specified. Model must be set by hand.");
   }
}

void ForeignInterfaceState::Evolve()
{
   if (_unloaded) {
      logger(Logger::ERROR, "Interface was previously unloaded. Cannot evolve.");
   }
   for (size_t i = 0; i < state.interm_scales.size(); i++) {
      state._solutions.emplace_back(Solution(state.discr, state._models, state.nf_initial_scale));

      // i=0: initial solution
      // i=1: O[0] S[0] = O(Q1 <- Q0) S(Q0)
      // i=2: O[1] S[1] = O[1] O[0] S[0] = O(Q2 <- Q1) O(Q1 <- Q0) S(Q0)
      for (int j = 0; j < static_cast<int>(i); j++) {
         ApplyEvolutionOperator(state._solutions[i], state.evol_op[j]);
      }
      state._solutions[i].RotateToPhysicalBasis();
   }

   for (const Solution &s : state._solutions) {
      state._fin_models.emplace_back(OutputModel(s));
   }
}

double ForeignInterfaceState::GetMoment(int which, double Q2, bool d2_or_elt)
{
   if (d2_or_elt && state.d2_weights == nullptr) {
      logger(Logger::WARNING, "GetMoment: trying to obtain `d2`, but the weights has not been "
                              "requested in the config. Ignoring request.");
      return NAN;
   }
   if (!d2_or_elt && state.elt_weights == nullptr) {
      logger(Logger::WARNING, "GetMoment: trying to obtain `ELT`, but the weights has not been "
                              "requested in the config. Ignoring request.");
      return NAN;
   }

   if (which < 1 || which > 6) {
      logger(Logger::ERROR, std::format("GetMoment: `which` should be the flavor index, "
                                        "and as such must be 1 <= which <= 6. I got: {:d}",
                                        which));
   }
   size_t flavor = static_cast<size_t>(which);
   size_t j      = 0;
   for (size_t i = 0; i < state.interm_scales.size(); i++, j++) {
      if (is_near(Q2, state.interm_scales[i], 1.0e-12)) {
         break;
      }
   }
   if (j == state.interm_scales.size()) {
      logger(Logger::ERROR, std::format("Scale: {:12e} is not among the available scales: ", Q2)
                                + vec_to_string(state.interm_scales));
   }
   if (flavor > state._solutions[j].nf) {
      logger(Logger::WARNING, "Requested flavor is not present at the requested scale...");
      return 0;
   }
   if (d2_or_elt) return state.d2_weights->ComputeSingleQuark(state._solutions[j]._distr_p[flavor]);
   else return state.elt_weights->ComputeSingleQuark(state._solutions[j]._distr_p[flavor]);
}

double ForeignInterfaceState::GetG2(double xBj, int which, double Q2)
{
   if (v_g2_weights == nullptr) {
      logger(Logger::WARNING, "GetG2: trying to obtain result, but the weights has not been "
                              "requested in the config. Ignoring request.");
      return NAN;
   }
   if (xBj < grid.grid_radius._d_info.intervals_phys[0].first || xBj > 1) {
      logger(Logger::WARNING, std::format("GetG2: trying to obtain result for xBj outside the "
                                          "range [{:.6f}, 1.0]. Ignoring request.",
                                          grid.grid_radius._d_info.intervals_phys[0].first));
      return NAN;
   }

   if (which < 0 || which > 6) {
      logger(Logger::ERROR, std::format("GetMoment: `which` should be the flavor index, "
                                        "and as such must be 1 <= which <= 6. I got: {:d}",
                                        which));
   }
   size_t flavor = static_cast<size_t>(which);
   size_t j      = 0;
   for (size_t i = 0; i < state.interm_scales.size(); i++, j++) {
      if (is_near(Q2, state.interm_scales[i], 1.0e-12)) {
         break;
      }
   }
   if (j == state.interm_scales.size()) {
      logger(Logger::ERROR, std::format("Scale: {:12e} is not among the available scales: ", Q2)
                                + vec_to_string(state.interm_scales));
   }
   if (flavor > state._solutions[j].nf) {
      logger(Logger::WARNING, "Requested flavor is not present at the requested scale...");
      return 0;
   }
   auto tmp = state.v_g2_weights->interpolate(xBj);

   return (tmp).dot(state._solutions[j]._distr_p[flavor]);
}

double ForeignInterfaceState::GetDistribution(OutputModel::FNC f, double Q2, double x1, double x2,
                                              double x3)
{
   size_t j = 0;
   for (size_t i = 0; i < state.interm_scales.size(); i++, j++) {
      if (is_near(Q2, state.interm_scales[i], 1.0e-12)) {
         break;
      }
   }
   if (j == state.interm_scales.size()) {
      logger(Logger::ERROR, std::format("Scale: {:12e} is not among the available scales: ", Q2)
                                + vec_to_string(state.interm_scales));
   }
   return state._fin_models[j].GetDistribution(f, x1, x2, x3);
}

void ForeignInterfaceState::Unload()
{
   state.evol_op.clear();
   state.evol_op.shrink_to_fit();
   state._unloaded = true;
}

}; // namespace Honeycomb

extern "C" {
void hc2_fi_set_up_(const char *config_name, int len)
{

   std::string conf_n(config_name, len);
   conf_n.erase(conf_n.find_last_not_of(' ') + 1);
   Honeycomb::logger(Honeycomb::Logger::INFO, "Loading config: \"" + conf_n + "\"");

   Honeycomb::timer::mark begin = Honeycomb::timer::now();
   Honeycomb::set_up(conf_n);
   Honeycomb::timer::mark end = Honeycomb::timer::now();
   Honeycomb::logger(
       Honeycomb::Logger::INFO,
       std::format("Full setup took {:.4e} (ms)", Honeycomb::timer::elapsed_ms(end, begin)));
}

void hc2_fi_set_model_(int *what, ModelSign model)
{
   if (*what < 0 || *what > 13) {
      Honeycomb::logger(Honeycomb::Logger::ERROR,
                        std::format("Model index {:d} is out of bound"
                                    " of available functions to set: [0, 13]",
                                    *what));
   }
   Honeycomb::InputModel::FNC f = static_cast<Honeycomb::InputModel::FNC>(*what);
   state._models.SetModel(f, [model](double x1, double x2, double x3) -> double {
      double a = x1, b = x2, c = x3;
      return model(&a, &b, &c);
   });
}

void hc2_fi_set_pim_model_(char *what, int len)
{
   std::string conf_n(what, len);
   conf_n.erase(conf_n.find_last_not_of(' ') + 1);
   Honeycomb::set_model(conf_n);
}

void hc2_fi_evolve_()
{
   state.Evolve();
}

void hc2_fi_unload_()
{
   state.Unload();
}

double hc2_fi_get_d2_(int *what, double *Q2)
{
   return state.GetMoment(*what, *Q2, true);
}

double hc2_fi_get_g2_(double *xBj, int *what, double *Q2)
{
   return state.GetG2(*xBj, *what, *Q2);
}

double hc2_fi_get_elt_(int *what, double *Q2)
{
   return state.GetMoment(*what, *Q2, false);
}

double hc2_fi_get_model_(int *what, double *Q2, double *x1, double *x2, double *x3)
{
   if (*what < 0 || *what > 13) {
      Honeycomb::logger(Honeycomb::Logger::ERROR,
                        std::format("Model index {:d} is out of bound"
                                    " of available functions to set: [0, 13]",
                                    *what));
   }
   Honeycomb::OutputModel::FNC f = static_cast<Honeycomb::OutputModel::FNC>(*what);
   return state.GetDistribution(f, *Q2, *x1, *x2, *x3);
}
}