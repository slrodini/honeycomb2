#include <honeycomb2/honeycomb2.hpp>
#include "honeycomb2/solution.hpp"
#include "honeycomb_120_25_grid_points.hpp"
#include <cereal/archives/json.hpp>

#include "original_model_functions.hpp"

double T_test(double x1, double x2, double x3)
{
   return (1 - x1 * x1) * (1 - x2 * x2) * (1 - x3 * x3);
}

double T_test_Hu(double x1, double x2, double x3)
{
   double r = std::max(std::fabs(x1), std::max(std::fabs(x2), std::fabs(x3)));
   return sin(x2 * M_PI) * (1 - x1 * x1) * (1 - x2 * x2) * (1 - x3 * x3) / (sqrt(r));
}

void test_save_and_laod();
void test_convolution();
void evolve_chiral_odd();
void evolve_chiral_even();
void test_solution_rotations();

int main()
{
   static_assert(Honeycomb::runge_kutta::RungeKuttaCompatible_2<Honeycomb::Kernels, Honeycomb::Solution>);

   // evolve_chiral_odd();
   // test_convolution();
   // test_save_and_laod();
   evolve_chiral_even();
   // test_solution_rotations();

   return 0;
}

void test_solution_rotations()
{

   const size_t n         = 6;
   const double rmin      = 0.01;
   Honeycomb::Grid2D grid = Honeycomb::generate_compliant_Grid2D(n, {rmin, 0.15, 0.65, 1}, {9, 9, 7});

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

   for (size_t nf = 1; nf <= 6; nf++) {
      Honeycomb::Solution sol0(&discr, model, nf);
      Honeycomb::Solution sol1(&discr, model, nf);
      sol1.RotateToPhysicalBasis();
      sol1.RotateToEvolutionBasis();

      double m = 0;
      for (size_t i = 0; i < sol1._distr_p.size(); i++) {
         for (long int k = 0; k < sol1._distr_p[i].size(); k++) {
            m = std::max(m, fabs(sol0._distr_p[i](k) - sol1._distr_p[i](k)));
            m = std::max(m, fabs(sol0._distr_m[i](k) - sol1._distr_m[i](k)));
         }
      }
      std::cout << std::format("{:.16e}\n", m);
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

   // auto [inter_scale, sol1]
   //     = get_initial_solution(Q02, Qf2, {0, 0, 0, 1.6129, 17.4724, 1.0e+6}, &discr, model);
   auto [inter_scale, sol1]
       = get_initial_solution(Q02, Qf2, {0, 0, 0, 1.0e+5, 2.0e+5, 1.0e+6}, &discr, model);

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

   auto callback = [&inter_scale](double t, const Honeycomb::Kernels &, Honeycomb::Solution &S) {
      if (t == inter_scale.back()) return;
      S.PushFlavor();
      return;
   };

   const double dt = 0.01;

   Honeycomb::runge_kutta::GenericRungeKutta<Honeycomb::Kernels, Honeycomb::Solution, 13> evolver(
       kers, sol1, Honeycomb::runge_kutta::DOPRI8, as, pref, t0, dt, callback);

   std::vector<double> back_scales;
   for (size_t i = inter_scale.size() - 2; i < inter_scale.size(); i--) {
      back_scales.push_back(inter_scale[i]);
   }
   back_scales.push_back(t0);

   fnc_elapsed = 0;
   begin       = Honeycomb::timer::now();
   // evolver(inter_scale, {10, 10, 40});
   evolver(inter_scale, std::vector<size_t>{40});
   end         = Honeycomb::timer::now();
   fnc_elapsed = Honeycomb::timer::elapsed_ms(end, begin);

   Honeycomb::logger(Honeycomb::Logger::INFO, std::format("  Elapsed: {:.4e} (ms)", fnc_elapsed));

   auto sol_fin = evolver.GetSolution();

   // auto callback_rev = [&back_scales](double t, const Honeycomb::Kernels &, Honeycomb::Solution &S) {
   //    if (t == back_scales.back()) return;
   //    S.PopFlavor();
   //    return;
   // };

   // Honeycomb::runge_kutta::GenericRungeKutta<Honeycomb::Kernels, Honeycomb::Solution, 13> evolver_back(
   //     kers, sol_fin, Honeycomb::runge_kutta::DOPRI8, as, pref, log(Qf2), dt, callback_rev);
   // evolver_back(back_scales, 20);

   // auto sol_bnf = evolver_back.GetSolution();

   std::vector<Honeycomb::Solution *> sols{&sol1, &sol_fin};

   for (auto s : sols) {
      double m_symm = 0;

      for (Honeycomb::RnC::Triplet &x123 : honeycomb_points) {

         Honeycomb::RnC::Pair rhophi_1_check = Honeycomb::RnC::from_x123_to_rhophi(x123);
         Honeycomb::RnC::Pair rhophi_2_check
             = Honeycomb::RnC::from_x123_to_rhophi(-x123[0], -x123[1], -x123[2]);

         // m_symm = std::max(m_symm, discr.interpolate_as_weights_v3(rhophi_1_check, s->_distr_p[0])
         //                               + discr.interpolate_as_weights_v3(rhophi_2_check, s->_distr_p[0]));

         // for (size_t k = 2; k < s->_distr_p.size(); k++) {
         //    m_symm = std::max(m_symm, discr.interpolate_as_weights_v3(rhophi_1_check, s->_distr_p[k])
         //                                  - discr.interpolate_as_weights_v3(rhophi_2_check,
         //                                  s->_distr_p[k]));
         // }

         // m_symm = std::max(m_symm, discr.interpolate_as_weights_v3(rhophi_1_check, s->_distr_m[0])
         //                               - discr.interpolate_as_weights_v3(rhophi_2_check,
         //                               s->_distr_m[0]));

         for (size_t k = 2; k < s->_distr_p.size(); k++) {
            m_symm = std::max(m_symm, discr.interpolate_as_weights_v3(rhophi_1_check, s->_distr_m[k])
                                          + discr.interpolate_as_weights_v3(rhophi_2_check, s->_distr_m[k]));
         }
      }
      Honeycomb::logger(Honeycomb::Logger::WARNING, std::format("Max symmetry violation: {:.16e}", m_symm));
   }

   std::FILE *fp_3 = std::fopen("CE_evo_check_evo_basis.dat", "w");

   for (Honeycomb::RnC::Triplet &x123 : honeycomb_points) {
      Honeycomb::RnC::Pair rhophi = Honeycomb::RnC::from_x123_to_rhophi(x123);

      std::fprintf(fp_3, "%.16e\t%.16e\t%.16e\t", x123[0], x123[1], x123[2]);
      for (size_t k = 0; k < sol_fin._distr_p.size(); k++) {
         std::fprintf(fp_3, "%.16e\t%.16e\t", discr.interpolate_as_weights_v3(rhophi, sol_fin._distr_p[k]),
                      discr.interpolate_as_weights_v3(rhophi, sol_fin._distr_m[k]));
      }
      std::fprintf(fp_3, "\n");
   }
   std::fclose(fp_3);

   sol1.RotateToPhysicalBasis();
   // sol_bnf.RotateToPhysicalBasis();
   sol_fin.RotateToPhysicalBasis();

   std::FILE *fp_2 = std::fopen("CE_evo_check.dat", "w");

   for (Honeycomb::RnC::Triplet &x123 : honeycomb_points) {
      Honeycomb::RnC::Pair rhophi = Honeycomb::RnC::from_x123_to_rhophi(x123);

      std::fprintf(fp_2, "%.16e\t%.16e\t%.16e\t", x123.v[0], x123.v[1], x123.v[2]);
      for (size_t k = 0; k < sol_fin._distr_p.size(); k++) {
         std::fprintf(fp_2, "%.16e\t%.16e\t", discr.interpolate_as_weights_v3(rhophi, sol_fin._distr_p[k]),
                      discr.interpolate_as_weights_v3(rhophi, sol_fin._distr_m[k]));
      }
      std::fprintf(fp_2, "\n");
   }
   std::fclose(fp_2);

   // double m_diff = 0;

   // for (Honeycomb::RnC::Triplet &x123 : honeycomb_points) {
   //    Honeycomb::RnC::Pair rhophi = Honeycomb::RnC::from_x123_to_rhophi(x123);
   //    for (size_t k = 0; k < sol_bnf._distr_p.size(); k++) {
   //       m_diff = std::max(m_diff, fabs(discr.interpolate_as_weights_v3(rhophi, sol_bnf._distr_p[k])
   //                                      - discr.interpolate_as_weights_v3(rhophi, sol1._distr_p[k])));
   //       m_diff = std::max(m_diff, fabs(discr.interpolate_as_weights_v3(rhophi, sol_bnf._distr_m[k])
   //                                      - discr.interpolate_as_weights_v3(rhophi, sol1._distr_m[k])));
   //    }
   // }
   // Honeycomb::logger(Honeycomb::Logger::INFO, std::format("Max difference: {:.10e}", m_diff));
}

void test_save_and_laod()
{
   const double Nc = 3;

   const size_t n         = 6;
   const double rmin      = 0.01;
   Honeycomb::Grid2D grid = Honeycomb::generate_compliant_Grid2D(n, {rmin, 0.15, 0.65, 1}, {9, 9, 7});

   // const size_t n         = 2;
   // const double rmin      = 0.01;
   // Honeycomb::Grid2D grid = Honeycomb::generate_compliant_Grid2D(n, {rmin, 0.15, 0.65, 1}, {2, 2, 2});

   double fnc_elapsed = 0;
   Honeycomb::logger(Honeycomb::Logger::INFO,
                     std::format("Total grid size: {:d}x{:d}", grid.size, grid.size));
   Honeycomb::logger(Honeycomb::Logger::INFO,
                     std::format("Total x123 size: {:d}x{:d}", grid.c_size, grid.c_size));

   // Discretize the kernels
   Honeycomb::timer::mark begin = Honeycomb::timer::now();
   Honeycomb::Kernels kers(grid, Nc);
   Honeycomb::timer::mark end = Honeycomb::timer::now();
   fnc_elapsed                = Honeycomb::timer::elapsed_ns(end, begin);
   Honeycomb::logger(Honeycomb::Logger::INFO,
                     std::format("Discretization time: {:+.6e} ms.", fnc_elapsed * 1.0e-6));

   // Save the discretized kernels
   begin = Honeycomb::timer::now();
   Honeycomb::save_kernels(kers, "check.cereal");

   end         = Honeycomb::timer::now();
   fnc_elapsed = Honeycomb::timer::elapsed_ns(end, begin);
   Honeycomb::logger(Honeycomb::Logger::INFO, std::format("Saving time: {:+.6e} ms.", fnc_elapsed * 1.0e-6));

   // Load the kernels
   begin                          = Honeycomb::timer::now();
   Honeycomb::Kernels kers_loaded = Honeycomb::load_kernels("check.cereal", grid, Nc);

   end         = Honeycomb::timer::now();
   fnc_elapsed = Honeycomb::timer::elapsed_ns(end, begin);
   Honeycomb::logger(Honeycomb::Logger::INFO, std::format("Loading time: {:+.6e} ms.", fnc_elapsed * 1.0e-6));

   for (long int i = 0; i < kers.H_NS.rows(); i++) {
      for (long int j = 0; j < kers.H_NS.cols(); j++) {
         if (std::fabs(kers.H_NS(i, j) - kers_loaded.H_NS(i, j)) > 2.0e-16) {
            Honeycomb::logger(Honeycomb::Logger::ERROR, "H_NS is not save/loaded correctly.");
         }
         if (std::fabs(kers.H_d13(i, j) - kers_loaded.H_d13(i, j)) > 2.0e-16) {
            Honeycomb::logger(Honeycomb::Logger::ERROR, "H_d13 is not save/loaded correctly.");
         }
         if (std::fabs(kers.H_gg_p(i, j) - kers_loaded.H_gg_p(i, j)) > 2.0e-16) {
            Honeycomb::logger(Honeycomb::Logger::ERROR, "H_gg_p is not save/loaded correctly.");
         }
         if (std::fabs(kers.H_gg_m(i, j) - kers_loaded.H_gg_m(i, j)) > 2.0e-16) {
            Honeycomb::logger(Honeycomb::Logger::ERROR, "H_gg_m is not save/loaded correctly.");
         }
         if (std::fabs(kers.H_qg_p(i, j) - kers_loaded.H_qg_p(i, j)) > 2.0e-16) {
            Honeycomb::logger(Honeycomb::Logger::ERROR, "H_qg_p is not save/loaded correctly.");
         }
         if (std::fabs(kers.H_qg_m(i, j) - kers_loaded.H_qg_m(i, j)) > 2.0e-16) {
            Honeycomb::logger(Honeycomb::Logger::ERROR, "H_qg_m is not save/loaded correctly.");
         }
         if (std::fabs(kers.H_gq_p(i, j) - kers_loaded.H_gq_p(i, j)) > 2.0e-16) {
            Honeycomb::logger(Honeycomb::Logger::ERROR, "H_gq_p is not save/loaded correctly.");
         }
         if (std::fabs(kers.H_gq_m(i, j) - kers_loaded.H_gq_m(i, j)) > 2.0e-16) {
            Honeycomb::logger(Honeycomb::Logger::ERROR, "H_gq_m is not save/loaded correctly.");
         }
      }
   }
}

void test_convolution()
{
   const double Nc = 3; // NC = 1 for tests

   const size_t n         = 6;
   const double rmin      = 0.01;
   Honeycomb::Grid2D grid = Honeycomb::generate_compliant_Grid2D(n, {rmin, 0.15, 0.65, 1}, {9, 9, 7});

   // const size_t n         = 13;
   // const double rmin      = 0.01;
   // Honeycomb::Grid2D grid = Honeycomb::generate_compliant_Grid2D(n, {rmin, 0.15, 0.65, 1}, {13, 13,
   // 11});

   double fnc_elapsed = 0;
   Honeycomb::logger(Honeycomb::Logger::INFO,
                     std::format("Total grid size: {:d}x{:d}", grid.size, grid.size));
   Honeycomb::logger(Honeycomb::Logger::INFO,
                     std::format("Total x123 size: {:d}x{:d}", grid.c_size, grid.c_size));

   Honeycomb::timer::mark begin = Honeycomb::timer::now();
   Honeycomb::Kernels kers(grid, Nc);
   Eigen::MatrixXd H_CO        = Honeycomb::get_CO_kernel(grid, Nc);
   Honeycomb::timer::mark end  = Honeycomb::timer::now();
   fnc_elapsed                += Honeycomb::timer::elapsed_ns(end, begin);

   Honeycomb::logger(Honeycomb::Logger::INFO,
                     std::format("Discretization time: {:+.6e} ms.", fnc_elapsed * 1.0e-6));

   const auto &KK = kers.H_NS;
   // const auto &KK = kers.H_test_1;

   for (long int i = 0; i < KK.rows(); i++) {
      for (long int j = 0; j < KK.cols(); j++) {
         if (std::isnan(KK(i, j)) || std::isinf(KK(i, j)))
            Honeycomb::logger(Honeycomb::Logger::WARNING,
                              std::format("NAN or INF in kernels! {:d}, {:d}", i, j));
      }
   }

   Honeycomb::Discretization discr(grid);
   Eigen::VectorXd fj    = discr(T_test);
   Eigen::VectorXd fj_Hu = discr(T_test_Hu);

   // Eigen::VectorXd fj2 = KK * fj;
   Eigen::VectorXd fj2 = H_CO * fj;

   std::FILE *fp = std::fopen("ConvCheck.dat", "w");
   for (long int i = 0; i < grid.c_size_li; i++) {
      std::fprintf(fp, "%.16e\t%.16e\t%.16e\t%.16e\t%.16e\n", grid._x123[i].v[0], grid._x123[i].v[1],
                   grid._x123[i].v[2], fj(i), fj2(i));
   }
   std::fclose(fp);

   double m1 = 0, m2 = 0;

   for (long int i = 0; i < grid.c_size_li; i++) {
      // -x123
      // auto rhophi = Honeycomb::RnC::from_x123_to_rhophi(
      //     Honeycomb::RnC::Triplet(-grid._x123[i].v[0], -grid._x123[i].v[1], -grid._x123[i].v[2]));

      // const double diff_1 = fj(i) - discr.interpolate_as_weights_v3(rhophi, fj);
      // const double diff_2 = fj2(i) + discr.interpolate_as_weights_v3(rhophi, fj2);
      auto rhophi = Honeycomb::RnC::from_x123_to_rhophi(
          Honeycomb::RnC::Triplet(-grid._x123[i].v[2], -grid._x123[i].v[1], -grid._x123[i].v[0]));

      const double diff_1 = fj(i) - discr.interpolate_as_weights_v3(rhophi, fj);
      const double diff_2 = fj2(i) - discr.interpolate_as_weights_v3(rhophi, fj2);

      m1 = std::max(m1, std::fabs(diff_1));
      m2 = std::max(m2, std::fabs(diff_2));
   }
   std::fprintf(stderr, "%.16e\t%.16e\n", m1, m2);
}

void evolve_chiral_odd()
{
   const double Nc = 3;

   const size_t n         = 6;
   const double rmin      = 0.01;
   Honeycomb::Grid2D grid = Honeycomb::generate_compliant_Grid2D(n, {rmin, 0.15, 0.65, 1}, {9, 9, 7});

   // const size_t n         = 13;
   // const double rmin      = 0.01;
   // Honeycomb::Grid2D grid = Honeycomb::generate_compliant_Grid2D(n, {rmin, 0.15, 0.65, 1}, {13, 13,
   // 11});

   double fnc_elapsed = 0;
   Honeycomb::logger(Honeycomb::Logger::INFO,
                     std::format("Total grid size: {:d}x{:d}", grid.size, grid.size));

   Honeycomb::timer::mark begin  = Honeycomb::timer::now();
   Eigen::MatrixXd H_CO          = Honeycomb::get_CO_kernel(grid, Nc);
   Honeycomb::timer::mark end    = Honeycomb::timer::now();
   fnc_elapsed                  += Honeycomb::timer::elapsed_ns(end, begin);

   Honeycomb::logger(Honeycomb::Logger::INFO,
                     std::format("Discretization time: {:+.6e} ms.", fnc_elapsed * 1.0e-6));

   Honeycomb::Discretization discr(grid);
   Eigen::VectorXd fj    = discr(T_test);
   Eigen::VectorXd fj_Hu = discr(T_test_Hu);

   const double Q02 = 1.0;
   const double Qf2 = 10000.0;

   const double t0 = log(Q02);
   const double tf = log(Qf2);

   const double dt_2 = 0.01; // ! Unused !

   const double pref = -1.0; // Evolution Equations are df/dt = - as Hxf

   Honeycomb::runge_kutta::GenericRungeKutta<Eigen::MatrixXd, Eigen::VectorXd, 13> evolver_CO(
       H_CO, fj, Honeycomb::runge_kutta::DOPRI8,
       [](double t) {
          return 1.0 / (11.0 * (t + 3.0));
       },
       pref, t0, dt_2,
       [](double, const Eigen::MatrixXd &, Eigen::VectorXd &) {
          return;
       });

   fnc_elapsed = 0;
   begin       = Honeycomb::timer::now();
   evolver_CO({tf}, 40);
   end          = Honeycomb::timer::now();
   fnc_elapsed += Honeycomb::timer::elapsed_ns(end, begin);

   auto fj_fin = evolver_CO.GetSolution();

   Honeycomb::logger(Honeycomb::Logger::INFO,
                     std::format("Evolution time: {:+.6e} ms.", fnc_elapsed * 1.0e-6));

   fnc_elapsed = 0;
   begin       = Honeycomb::timer::now();
   evolver_CO.reset(fj_Hu, t0);
   evolver_CO({log(10000)}, 40);
   end          = Honeycomb::timer::now();
   fnc_elapsed += Honeycomb::timer::elapsed_ns(end, begin);

   Honeycomb::logger(Honeycomb::Logger::INFO,
                     std::format("Evolution time: {:+.6e} ms.", fnc_elapsed * 1.0e-6));

   auto fj_Hu_fin = evolver_CO.GetSolution();

   std::FILE *fp_2 = std::fopen("CO_evo_check.dat", "w");

   for (Honeycomb::RnC::Triplet &x123 : honeycomb_points) {
      Honeycomb::RnC::Pair rhophi = Honeycomb::RnC::from_x123_to_rhophi(x123);
      std::fprintf(
          fp_2, "%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\n", x123.v[0], x123.v[1], x123.v[2],
          discr.interpolate_as_weights_v3(rhophi, fj), discr.interpolate_as_weights_v3(rhophi, fj_fin),
          discr.interpolate_as_weights_v3(rhophi, fj_Hu), discr.interpolate_as_weights_v3(rhophi, fj_Hu_fin));
   }
   std::fclose(fp_2);
}