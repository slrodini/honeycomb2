#include <honeycomb2/honeycomb2.hpp>
#include "honeycomb2/solution.hpp"
#include "honeycomb_120_25_grid_points.hpp"
#include <cereal/archives/json.hpp>

#include "original_model_functions.hpp"

void test_convolution();
void test_rk_step();
void evolve_chiral_even();

int main()
{
   static_assert(Honeycomb::runge_kutta::RungeKuttaCompatible_2<Honeycomb::Kernels, Honeycomb::Solution>);

   // test_convolution();
   evolve_chiral_even();

   return 0;
}

void test_convolution()
{

   Honeycomb::timer::mark begin = Honeycomb::timer::now();
   Honeycomb::timer::mark end   = Honeycomb::timer::now();

   // fnc_elapsed                = Honeycomb::timer::elapsed_ns(end, begin);

   const double Nc = 3; // NC = 1 for tests

   const size_t n         = 8;
   const double rmin      = 0.01;
   Honeycomb::Grid2D grid = Honeycomb::generate_compliant_Grid2D(n, {rmin, 0.1, 0.4, 1}, {12, 9, 7});

   Honeycomb::logger(Honeycomb::Logger::INFO,
                     std::format("Total grid size: {:d}x{:d}", grid.size, grid.size));
   Honeycomb::logger(Honeycomb::Logger::INFO,
                     std::format("Total x123 size: {:d}x{:d}", grid.c_size, grid.c_size));

   Honeycomb::logger(Honeycomb::Logger::INFO, std::format("Kernel Discretization"));
   begin = Honeycomb::timer::now();
   Honeycomb::Kernels kers(grid, Nc);
   // Honeycomb::save_kernels<cereal::PortableBinaryOutputArchive>(kers, "n6Ker.cereal");
   // Honeycomb::Kernels kers = Honeycomb::load_kernels("kers_for_checks.cereal", grid, Nc);

   auto &KK = kers.H_qg_p;
   for (long int i = 0; i < KK.rows(); i++) {
      for (long int j = 0; j < KK.cols(); j++) {
         if (std::isnan(KK(i, j)) || std::isinf(KK(i, j)))
            Honeycomb::logger(Honeycomb::Logger::WARNING,
                              std::format("NAN or INF in kernels! {:d}, {:d}", i, j));
      }
   }

   end = Honeycomb::timer::now();
   Honeycomb::logger(Honeycomb::Logger::INFO,
                     std::format("  Elapsed: {:.4e} (ms)", Honeycomb::timer::elapsed_ms(end, begin)));

   Honeycomb::Discretization discr(grid);

   Honeycomb::InputModel model;

   auto fnc_sym = [](double x1, double x2, double x3) {
      return (1 - x1 * x1) * (1 - x2 * x2) * (1 - x3 * x3);
   };

   auto a_fnc_sym = [](double x1, double x2, double x3) {
      return x2 * (1 - x1 * x1) * (1 - x2 * x2) * (1 - x3 * x3);
   };

   model.SetModel(Honeycomb::InputModel::T_UP, fnc_sym);
   model.SetModel(Honeycomb::InputModel::DT_UP, a_fnc_sym);

   model.SetModel(Honeycomb::InputModel::T_DN, fnc_sym);
   model.SetModel(Honeycomb::InputModel::DT_DN, a_fnc_sym);

   model.SetModel(Honeycomb::InputModel::T_ST, fnc_sym);
   model.SetModel(Honeycomb::InputModel::DT_ST, a_fnc_sym);

   model.SetModel(Honeycomb::InputModel::T_P_GL, fnc_sym);
   model.SetModel(Honeycomb::InputModel::T_M_GL, a_fnc_sym);

   const double Q02 = 1.0;
   const double Qf2 = 1.0e+4;

   // auto [inter_scale, sol1]
   //     = get_initial_solution(Q02, Qf2, {0, 0, 0, 1.6129, 17.4724, 1.0e+6}, &discr, model);
   auto [inter_scale, sol1]
       = get_initial_solution(Q02, Qf2, {0, 0, 0, 1.0e+5, 2.0e+5, 1.0e+6}, &discr, model);
   {
      std::FILE *fp_1 = std::fopen("Initial_Solution.dat", "w");

      for (Honeycomb::RnC::Triplet &x123 : honeycomb_points) {
         Honeycomb::RnC::Pair rhophi = Honeycomb::RnC::from_x123_to_rhophi(x123);

         std::fprintf(fp_1, "%.16e\t%.16e\t%.16e\t", x123[0], x123[1], x123[2]);
         for (size_t k = 0; k < sol1._distr_p.size(); k++) {
            std::fprintf(fp_1, "%.16e\t%.16e\t", discr.interpolate_as_weights_v3(rhophi, sol1._distr_p[k]),
                         discr.interpolate_as_weights_v3(rhophi, sol1._distr_m[k]));
         }
         std::fprintf(fp_1, "\n");
      }
      std::fclose(fp_1);
   }

   sol1._ker_mul(-1.0, kers);

   {
      std::FILE *fp_1 = std::fopen("Convolved_Solution.dat", "w");

      for (Honeycomb::RnC::Triplet &x123 : honeycomb_points) {
         Honeycomb::RnC::Pair rhophi = Honeycomb::RnC::from_x123_to_rhophi(x123);

         std::fprintf(fp_1, "%.16e\t%.16e\t%.16e\t", x123[0], x123[1], x123[2]);
         for (size_t k = 0; k < sol1._distr_p.size(); k++) {
            std::fprintf(fp_1, "%.16e\t%.16e\t", discr.interpolate_as_weights_v3(rhophi, sol1._distr_p[k]),
                         discr.interpolate_as_weights_v3(rhophi, sol1._distr_m[k]));
         }
         std::fprintf(fp_1, "\n");
      }
      std::fclose(fp_1);
   }
}

void test_rk_step()
{

   Honeycomb::timer::mark begin = Honeycomb::timer::now();
   Honeycomb::timer::mark end   = Honeycomb::timer::now();

   // fnc_elapsed                = Honeycomb::timer::elapsed_ns(end, begin);

   const double Nc = 3; // NC = 1 for tests

   const size_t n         = 8;
   const double rmin      = 0.01;
   Honeycomb::Grid2D grid = Honeycomb::generate_compliant_Grid2D(n, {rmin, 0.1, 0.4, 1}, {12, 9, 7});

   Honeycomb::logger(Honeycomb::Logger::INFO,
                     std::format("Total grid size: {:d}x{:d}", grid.size, grid.size));
   Honeycomb::logger(Honeycomb::Logger::INFO,
                     std::format("Total x123 size: {:d}x{:d}", grid.c_size, grid.c_size));

   Honeycomb::logger(Honeycomb::Logger::INFO, std::format("Kernel Discretization"));
   begin = Honeycomb::timer::now();
   // Honeycomb::Kernels kers(grid, Nc);
   // Honeycomb::save_kernels<cereal::PortableBinaryOutputArchive>(kers, "n6Ker.cereal");
   Honeycomb::Kernels kers = Honeycomb::load_kernels("kers_for_checks.cereal", grid, Nc);

   end = Honeycomb::timer::now();
   Honeycomb::logger(Honeycomb::Logger::INFO,
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

   const double Q02 = 1.0;
   const double Qf2 = 1.0e+4;

   // auto [inter_scale, sol1]
   //     = get_initial_solution(Q02, Qf2, {0, 0, 0, 1.6129, 17.4724, 1.0e+6}, &discr, model);
   auto [inter_scale, sol1]
       = get_initial_solution(Q02, Qf2, {0, 0, 0, 1.0e+5, 2.0e+5, 1.0e+6}, &discr, model);
   {
      std::FILE *fp_1 = std::fopen("Initial_Solution.dat", "w");

      for (Honeycomb::RnC::Triplet &x123 : honeycomb_points) {
         Honeycomb::RnC::Pair rhophi = Honeycomb::RnC::from_x123_to_rhophi(x123);

         std::fprintf(fp_1, "%.16e\t%.16e\t%.16e\t", x123[0], x123[1], x123[2]);
         for (size_t k = 0; k < sol1._distr_p.size(); k++) {
            std::fprintf(fp_1, "%.16e\t%.16e\t", discr.interpolate_as_weights_v3(rhophi, sol1._distr_p[k]),
                         discr.interpolate_as_weights_v3(rhophi, sol1._distr_m[k]));
         }
         std::fprintf(fp_1, "\n");
      }
      std::fclose(fp_1);
   }

   auto as = [](double t) -> double {
      return 1.0 / (11.0 * (t + 3.0));
   };
   const double pref = -1.0;
   const double t0   = log(Q02);
   const double dt   = 0.25;
   Honeycomb::runge_kutta::GenericRungeKutta<Honeycomb::Kernels, Honeycomb::Solution, 13> evolver(
       kers, sol1, Honeycomb::runge_kutta::DOPRI8, as, pref, t0, dt);

   evolver();

   auto sol2 = evolver.GetSolution();

   {
      std::FILE *fp_1 = std::fopen("Convolved_Solution.dat", "w");

      for (Honeycomb::RnC::Triplet &x123 : honeycomb_points) {
         Honeycomb::RnC::Pair rhophi = Honeycomb::RnC::from_x123_to_rhophi(x123);

         std::fprintf(fp_1, "%.16e\t%.16e\t%.16e\t", x123[0], x123[1], x123[2]);
         for (size_t k = 0; k < sol2._distr_p.size(); k++) {
            std::fprintf(fp_1, "%.16e\t%.16e\t", discr.interpolate_as_weights_v3(rhophi, sol2._distr_p[k]),
                         discr.interpolate_as_weights_v3(rhophi, sol2._distr_m[k]));
         }
         std::fprintf(fp_1, "\n");
      }
      std::fclose(fp_1);
   }
}

void evolve_chiral_even()
{

   double fnc_elapsed           = 0;
   Honeycomb::timer::mark begin = Honeycomb::timer::now();
   Honeycomb::timer::mark end   = Honeycomb::timer::now();

   // fnc_elapsed                = Honeycomb::timer::elapsed_ns(end, begin);

   const double Nc = 3; // NC = 1 for tests

   const size_t n         = 8;
   const double rmin      = 0.01;
   Honeycomb::Grid2D grid = Honeycomb::generate_compliant_Grid2D(n, {rmin, 0.1, 0.4, 1}, {12, 9, 7});

   Honeycomb::logger(Honeycomb::Logger::INFO,
                     std::format("Total grid size: {:d}x{:d}", grid.size, grid.size));
   Honeycomb::logger(Honeycomb::Logger::INFO,
                     std::format("Total x123 size: {:d}x{:d}", grid.c_size, grid.c_size));

   Honeycomb::logger(Honeycomb::Logger::INFO, std::format("Kernel Discretization"));
   begin = Honeycomb::timer::now();
   // Honeycomb::Kernels kers(grid, Nc);
   // Honeycomb::save_kernels<cereal::PortableBinaryOutputArchive>(kers, "n6Ker.cereal");
   Honeycomb::Kernels kers = Honeycomb::load_kernels("n6Ker.cereal", grid, Nc);
   Honeycomb::MergedKernelsFixedNf merge_ker(kers, 3);

   end = Honeycomb::timer::now();
   Honeycomb::logger(Honeycomb::Logger::INFO,
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
   const double Qf2 = 1.0e+4;

   const double t0 = log(Q02);
   const double tf = log(Qf2);

   Honeycomb::Solution sol1(&discr, model, 3);
   Honeycomb::EvolutionOperatorFixedNf O1(&grid, 3);

   {
      std::FILE *fp_1 = std::fopen("CE_evo_check_evo_basis_initial.dat", "w");

      for (Honeycomb::RnC::Triplet &x123 : honeycomb_points) {
         Honeycomb::RnC::Pair rhophi = Honeycomb::RnC::from_x123_to_rhophi(x123);

         std::fprintf(fp_1, "%.16e\t%.16e\t%.16e\t", x123[0], x123[1], x123[2]);
         for (size_t k = 0; k < sol1._distr_p.size(); k++) {
            std::fprintf(fp_1, "%.16e\t%.16e\t", discr.interpolate_as_weights_v3(rhophi, sol1._distr_p[k]),
                         discr.interpolate_as_weights_v3(rhophi, sol1._distr_m[k]));
         }
         std::fprintf(fp_1, "\n");
      }
      std::fclose(fp_1);
   }

   const double dt = 0.01;

   Honeycomb::runge_kutta::GenericRungeKutta<Honeycomb::Kernels, Honeycomb::Solution, 13> evolver(
       kers, sol1, Honeycomb::runge_kutta::DOPRI8, as, pref, t0, dt);

   Honeycomb::runge_kutta::GenericRungeKutta<Honeycomb::MergedKernelsFixedNf,
                                             Honeycomb::EvolutionOperatorFixedNf, 13>
       evolver_O(merge_ker, O1, Honeycomb::runge_kutta::DOPRI8, as, pref, t0, dt);

   fnc_elapsed = 0;
   begin       = Honeycomb::timer::now();
   evolver({tf}, std::vector<size_t>{40});
   end         = Honeycomb::timer::now();
   fnc_elapsed = Honeycomb::timer::elapsed_ms(end, begin);

   Honeycomb::logger(Honeycomb::Logger::INFO, std::format("  Elapsed: {:.4e} (ms)", fnc_elapsed));

   begin = Honeycomb::timer::now();
   evolver_O({tf}, std::vector<size_t>{40});
   // evolver_O();
   end         = Honeycomb::timer::now();
   fnc_elapsed = Honeycomb::timer::elapsed_ms(end, begin);
   Honeycomb::logger(Honeycomb::Logger::INFO, std::format("  Elapsed: {:.4e} (ms)", fnc_elapsed));

   auto sol_fin = evolver.GetSolution();
   auto O_fin   = evolver_O.GetSolution();
   auto sol2    = sol1;
   Honeycomb::ApplyEvolutionOperator(sol2, O_fin);

   std::vector<Honeycomb::Solution *> sols{&sol1, &sol_fin, &sol2};

   for (auto s : sols) {
      double m_symm = 0;

      for (Honeycomb::RnC::Triplet &x123 : honeycomb_points) {

         Honeycomb::RnC::Pair rhophi_1_check = Honeycomb::RnC::from_x123_to_rhophi(x123);
         Honeycomb::RnC::Pair rhophi_2_check
             = Honeycomb::RnC::from_x123_to_rhophi(-x123[0], -x123[1], -x123[2]);

         for (size_t k = 2; k < s->_distr_p.size(); k++) {
            m_symm = std::max(m_symm, discr.interpolate_as_weights_v3(rhophi_1_check, s->_distr_m[k])
                                          + discr.interpolate_as_weights_v3(rhophi_2_check, s->_distr_m[k]));
         }
      }
      Honeycomb::logger(Honeycomb::Logger::WARNING, std::format("Max symmetry violation: {:.16e}", m_symm));
   }

   // sol1.RotateToPhysicalBasis();
   // sol_fin.RotateToPhysicalBasis();
   // sol2.RotateToPhysicalBasis();

   std::FILE *fp_2 = std::fopen("CE_evo_check.dat", "w");

   for (Honeycomb::RnC::Triplet &x123 : honeycomb_points) {
      Honeycomb::RnC::Pair rhophi = Honeycomb::RnC::from_x123_to_rhophi(x123);

      std::fprintf(fp_2, "%.16e\t%.16e\t%.16e\t", x123.v[0], x123.v[1], x123.v[2]);
      for (size_t k = 0; k < sol_fin._distr_p.size(); k++) {
         std::fprintf(fp_2, "%.16e\t%.16e\t%.16e\t%.16e\t",
                      discr.interpolate_as_weights_v3(rhophi, sol_fin._distr_p[k]),
                      discr.interpolate_as_weights_v3(rhophi, sol_fin._distr_m[k]),
                      discr.interpolate_as_weights_v3(rhophi, sol2._distr_p[k]),
                      discr.interpolate_as_weights_v3(rhophi, sol2._distr_m[k]));
      }
      std::fprintf(fp_2, "\n");
   }
   std::fclose(fp_2);
}
