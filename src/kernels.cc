#include <honeycomb2/kernels.hpp>
#include <honeycomb2/thread_pool.hpp>

namespace Honeycomb
{

void InputModel::SetModel(FNC f, std::function<double(double, double, double)> model)
{
   switch (f) {
   case T_UP:
      T_up = model;
      break;
   case T_DN:
      T_dn = model;
      break;
   case T_ST:
      T_st = model;
      break;
   case T_CH:
      T_ch = model;
      break;
   case T_BM:
      T_bm = model;
      break;
   case T_TP:
      T_tp = model;
      break;
   case DT_UP:
      DT_up = model;
      break;
   case DT_DN:
      DT_dn = model;
      break;
   case DT_ST:
      DT_st = model;
      break;
   case DT_CH:
      DT_ch = model;
      break;
   case DT_BM:
      DT_bm = model;
      break;
   case DT_TP:
      DT_tp = model;
      break;
   case T_P_GL:
      T_p_gl = model;
      break;
   case T_M_GL:
      T_m_gl = model;
      break;
   default:
      break;
   }
}

Eigen::MatrixXd get_CO_kernel(const Grid2D &g, double _Nc)
{
   const double CF = (_Nc * _Nc - 1.0) / (2.0 * _Nc);

   Eigen::MatrixXd H_CO = Eigen::MatrixXd::Zero(g.c_size, g.c_size);

   ThreadPool th_pool(10);
   // std::mutex th_mutex;

   for (size_t c_a = 0; c_a < g.c_size; c_a++) {
      // Unsubtracted integrals
      th_pool.AddTask([&, c_a](void) -> void {
         std::unique_lock<std::mutex> lock(th_pool.task_mutex, std::defer_lock);

         for (size_t aP = 0; aP < g.size; aP++) {
            auto [jP, iP] = g.get_double_index(aP);
            size_t c_iP   = g.grid_angle._from_iw_to_ic[iP];
            size_t c_jP   = g.grid_radius._from_iw_to_ic[jP];
            size_t c_aP   = g.c_get_flatten_index(c_jP, c_iP);

            double res_plus_12 = Hplus12::integrate(c_a, aP, g);
            double res_plus_23 = Hplus23::integrate(c_a, aP, g);

            double res_minus_12 = Hminus12::integrate(c_a, aP, g);
            double res_minus_23 = Hminus23::integrate(c_a, aP, g);

            lock.lock();

            H_CO(c_a, c_aP) -= _Nc * 2.0 * (res_plus_12 + res_plus_23);
            H_CO(c_a, c_aP) -= (2.0 / _Nc) * (res_minus_12 + res_minus_23);

            lock.unlock();
         }
      });

      // Subtracted integrals
      th_pool.AddTask([&, c_a](void) -> void {
         std::unique_lock<std::mutex> lock(th_pool.task_mutex, std::defer_lock);
         for (size_t c_aP = 0; c_aP < g.c_size; c_aP++) {
            double res_hat_12 = Hhat12::subtracted_integrate(c_a, c_aP, g);
            double res_hat_23 = Hhat23::subtracted_integrate(c_a, c_aP, g);
            double res_hat_13 = Hhat13::subtracted_integrate(c_a, c_aP, g);

            lock.lock();

            H_CO(c_a, c_aP) += _Nc * (res_hat_12 + res_hat_23);
            H_CO(c_a, c_aP) += (1.0 / _Nc) * (-res_hat_13);
            lock.unlock();
         }
      });
   }

   th_pool.WaitOnJobs();

   // Field anomalous dimensions, nf independent
   for (size_t c_a = 0; c_a < g.c_size; c_a++) {
      H_CO(c_a, c_a) -= 3.0 * CF;
   }

   return H_CO;
}

Kernels::Kernels(const Grid2D &g, double _Nc, bool to_compute)
    : grid(g), Nc(_Nc), CA(_Nc), CF((_Nc * _Nc - 1) / (2.0 * _Nc))
{

   // NOTE: The field anomalous dimensions ARE NOT included in the kernels.
   H_NS   = Eigen::MatrixXd::Zero(g.c_size, g.c_size);
   H_d13  = Eigen::MatrixXd::Zero(g.c_size, g.c_size);
   H_gg_p = Eigen::MatrixXd::Zero(g.c_size, g.c_size);
   H_gg_m = Eigen::MatrixXd::Zero(g.c_size, g.c_size);
   H_qg_p = Eigen::MatrixXd::Zero(g.c_size, g.c_size);
   H_qg_m = Eigen::MatrixXd::Zero(g.c_size, g.c_size);
   H_gq_p = Eigen::MatrixXd::Zero(g.c_size, g.c_size);
   H_gq_m = Eigen::MatrixXd::Zero(g.c_size, g.c_size);

   H_test_1 = Eigen::MatrixXd::Zero(g.c_size, g.c_size);
   H_test_2 = Eigen::MatrixXd::Zero(g.c_size, g.c_size);

   if (!to_compute) return;

   ThreadPool th_pool(10);
   // std::mutex th_mutex;

   for (size_t c_a = 0; c_a < g.c_size; c_a++) {
      // Unsubtracted integrals
      th_pool.AddTask([&, c_a](void) -> void {
         std::unique_lock<std::mutex> lock(th_pool.task_mutex, std::defer_lock);

         for (size_t aP = 0; aP < g.size; aP++) {
            auto [jP, iP] = g.get_double_index(aP);
            size_t c_iP   = g.grid_angle._from_iw_to_ic[iP];
            size_t c_jP   = g.grid_radius._from_iw_to_ic[jP];
            size_t c_aP   = g.c_get_flatten_index(c_jP, c_iP);

            // // NS
            double res_plus_12  = Hplus12::integrate(c_a, aP, g);
            double res_plus_13  = Hplus13::integrate(c_a, aP, g);
            double res_minus_12 = Hminus12::integrate(c_a, aP, g);
            double res_he23p23  = He23P23::integrate(c_a, aP, g);

            // qq S => NO nf here, multipliead in the evolution
            double res_d_13 = 4.0 * Hd13::integrate(c_a, aP, g);

            // GG
            // This is -4H^+ -2\tilde{H}^+
            double res_plus_12_gg = HplusSum12GG::integrate(c_a, aP, g);
            double res_plus_13_gg = HplusSum13GG::integrate(c_a, aP, g);

            double res_minus_12_gg = Hminus12GG::integrate(c_a, aP, g);
            double res_minus_13_gg = Hminus13GG::integrate(c_a, aP, g);

            // QG
            double res_v_minus = Vminus13::integrate(c_a, aP, g);

            double res_v_plus = Vplus13::integrate(c_a, aP, g); // For testing purposes, to be removed

            // GQ
            double res_w_plus     = Wplus13::integrate(c_a, aP, g);
            double res_w_plus_p23 = Wplus13P23::integrate(c_a, aP, g);

            double res_w_minus     = Wminus13::integrate(c_a, aP, g);
            double res_w_minus_p23 = Wminus13P23::integrate(c_a, aP, g);

            double res_dw     = DeltaW::integrate(c_a, aP, g);
            double res_dw_p23 = DeltaWP23::integrate(c_a, aP, g);

            lock.lock();
            H_NS(c_a, c_aP) += _Nc * 2.0 * (-res_plus_12);
            H_NS(c_a, c_aP) += (1.0 / _Nc) * (res_plus_13 - 2.0 * res_minus_12 + res_he23p23);

            H_d13(c_a, c_aP) += res_d_13;

            H_gg_p(c_a, c_aP) += _Nc * (res_plus_12_gg + res_plus_13_gg + res_minus_12_gg + res_minus_13_gg);
            H_gg_m(c_a, c_aP) += _Nc * (res_plus_12_gg + res_plus_13_gg - res_minus_12_gg - res_minus_13_gg);

            // Note: NO nf here, multiplied in the evolution
            H_qg_p(c_a, c_aP) += -res_v_minus;
            H_qg_m(c_a, c_aP) += +res_v_minus;

            H_test_2(c_a, c_aP) += res_v_plus;

            H_gq_p(c_a, c_aP) +=
                _Nc * (res_w_plus - res_w_plus_p23 + res_w_minus - res_w_minus_p23 - 2 * res_dw + 2 * res_dw_p23);
            H_gq_m(c_a, c_aP) += -(_Nc - 4.0 / _Nc) * (res_w_plus + res_w_plus_p23 + res_w_minus + res_w_minus_p23);

            lock.unlock();
         }
      });

      // Subtracted integrals
      th_pool.AddTask([&, c_a](void) -> void {
         std::unique_lock<std::mutex> lock(th_pool.task_mutex, std::defer_lock);
         for (size_t c_aP = 0; c_aP < g.c_size; c_aP++) {

            // NS
            double res_hat_12 = Hhat12::subtracted_integrate(c_a, c_aP, g);
            double res_hat_23 = Hhat23::subtracted_integrate(c_a, c_aP, g);
            double res_hat_13 = Hhat13::subtracted_integrate(c_a, c_aP, g);

            // GG
            double res_hat_12_gg = Hhat12GG::subtracted_integrate(c_a, c_aP, g);
            double res_hat_23_gg = Hhat23GG::subtracted_integrate(c_a, c_aP, g);
            double res_hat_31_gg = Hhat31GG::subtracted_integrate(c_a, c_aP, g);

            // QG
            double res_v_plus = Vplus13::subtracted_integrate(c_a, c_aP, g);

            lock.lock();
            H_NS(c_a, c_aP) += _Nc * (res_hat_12 + res_hat_23);
            H_NS(c_a, c_aP) += (1.0 / _Nc) * (-res_hat_13);

            H_gg_p(c_a, c_aP) += _Nc * (res_hat_12_gg + res_hat_23_gg + res_hat_31_gg);
            H_gg_m(c_a, c_aP) += _Nc * (res_hat_12_gg + res_hat_23_gg + res_hat_31_gg);

            H_qg_p(c_a, c_aP) += res_v_plus;
            H_qg_m(c_a, c_aP) += res_v_plus;

            H_test_1(c_a, c_aP) += res_v_plus;

            lock.unlock();
         }
      });
   }

   th_pool.WaitOnJobs();
};

} // namespace Honeycomb
