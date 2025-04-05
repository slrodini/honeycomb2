#include <honeycomb2/kernels.hpp>
#include <honeycomb2/thread_pool.hpp>

// TODO: The logic seems sound becuase Hplus12 convolved with a test function is prefectly smooth
// TODO: The problem is in Hhat12, which is likely the plus prescription that is giving issues.

namespace Honeycomb
{

Kernels::Kernels(const Grid2D &g, double _Nc) : grid(g), Nc(_Nc), CA(_Nc), CF((_Nc * _Nc - 1) / (2.0 * _Nc))
{

   ThreadPool th_pool(10);
   std::mutex th_mutex;
   std::condition_variable th_done;

   H_NS   = Eigen::MatrixXd::Zero(g.c_size, g.c_size);
   H_CO   = Eigen::MatrixXd::Zero(g.c_size, g.c_size);
   H_d13  = Eigen::MatrixXd::Zero(g.c_size, g.c_size);
   H_gg_p = Eigen::MatrixXd::Zero(g.c_size, g.c_size);
   H_gg_m = Eigen::MatrixXd::Zero(g.c_size, g.c_size);
   H_qg_p = Eigen::MatrixXd::Zero(g.c_size, g.c_size);
   H_qg_m = Eigen::MatrixXd::Zero(g.c_size, g.c_size);
   H_gq_p = Eigen::MatrixXd::Zero(g.c_size, g.c_size);
   H_gq_m = Eigen::MatrixXd::Zero(g.c_size, g.c_size);

   H_plus_12_v1 = Eigen::MatrixXd::Zero(g.c_size, g.c_size);
   H_plus_12_v2 = Eigen::MatrixXd::Zero(g.c_size, g.c_size);

   // H_{aa'} = H_{[c_i c_j],[c_i' c_j']} where c_i over angles and c_j over radii
   // Then c_. can be obtained from weight index via mapping,
   // So, I loop over weight index to have simple expressions for integrals

   // Timings: for the grid in twodim_check the following code takes about 28 seconds.
   // Maybe I can improve it using the thread_pool.
   // With thread_pool in inner loop around 44 seconds
   // With thread_pool on outer loop (i.e. all inner loop in one task), around 26 seconds

   size_t tot_task = g.c_size + g.c_size;
   size_t num_done = 0;

   for (size_t c_a = 0; c_a < g.c_size; c_a++) {
      // Unsubtracted integrals
      th_pool.add_task([&, c_a](void) -> void {
         std::unique_lock<std::mutex> lock(th_mutex, std::defer_lock);

         for (size_t aP = 0; aP < g.size; aP++) {
            auto [jP, iP] = g.get_double_index(aP);
            size_t c_iP   = g.grid_angle._from_iw_to_ic[iP];
            size_t c_jP   = g.grid_radius._from_iw_to_ic[jP];
            size_t c_aP   = g.c_get_flatten_index(c_jP, c_iP);

            double res = Hplus12::integrate(c_a, aP, g);

            lock.lock();
            H_plus_12_v1(c_a, c_aP) += res;
            H_NS(c_a, c_aP)         += res;

            lock.unlock();
         }
         lock.lock();
         num_done++;
         lock.unlock();
         th_done.notify_one();
      });

      // Subtracted integrals
      th_pool.add_task([&, c_a](void) -> void {
         std::unique_lock<std::mutex> lock(th_mutex, std::defer_lock);
         for (size_t c_aP = 0; c_aP < g.c_size; c_aP++) {
            double res  = Hhat12::subtracted_integrate(c_a, c_aP, g);
            double res2 = Hplus12::integrate_v2(c_a, c_aP, g);

            lock.lock();
            H_plus_12_v2(c_a, c_aP) += res2;
            H_NS(c_a, c_aP)         += res;
            lock.unlock();
         }
         lock.lock();
         num_done++;
         lock.unlock();
         th_done.notify_one();
      });
   }

   std::unique_lock<std::mutex> lock_done(th_mutex);
   th_done.wait(lock_done, [&num_done, &tot_task](void) {
      return num_done >= tot_task;
   });

   // Here compute Integrate H(x123_[c_i c_j], v) \times w_{ij}(\rho, \phi)
   // transformation of w_{ij} is done in constructor.
}; // namespace Honeycomb

} // namespace Honeycomb