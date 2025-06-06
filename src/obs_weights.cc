#include "gauss_kronrod.hpp"
#include <honeycomb2/obs_weights.hpp>
#include <honeycomb2/thread_pool.hpp>
#include <honeycomb2/utilities.hpp>

namespace Honeycomb
{

Eigen::VectorXd G2Weights::GetWeights(double xBj, const Grid2D &grid, double int_e_r, double int_e_a)
{

   if (!grid.is_compliant) {
      logger(Logger::ERROR, "G2Weights: The Grid2D is not compliant to the specs. Please use only grids "
                            "generated via `generate_compliant_Grid2D`.");
   }
   using integrator            = Honeycomb::Integration::GaussKronrod<Honeycomb::Integration::GK_61>;
   unsigned int num_av_threads = std::thread::hardware_concurrency();
   unsigned int num_threads    = num_av_threads <= 2 ? 1 : num_av_threads - 2;
   ThreadPool th(num_threads);

   std::mutex done_mutex;

   Eigen::VectorXd weights = Eigen::VectorXd::Zero(grid.c_size_li);

   //  Sector 0, only first integral
   for (size_t k = grid.grid_angle._delim_indexes[0]; k < grid.grid_angle._delim_indexes[1]; k++) {
      for (size_t j = 0; j < grid.grid_radius.size; j++) {
         th.AddTask([&, j, k](void) -> void {
            auto r_supp         = grid.grid_radius.get_support_weight_aj(j);
            const double xi_min = std::max(grid.grid_radius._d_info.to_phys_space(r_supp.first), xBj);
            const double xi_max = std::min(grid.grid_radius._d_info.to_phys_space(r_supp.second), 1.0);
            auto index          = grid.get_flatten_index(j, k);
            size_t c_k          = grid.grid_angle._from_iw_to_ic[k];
            size_t c_j          = grid.grid_radius._from_iw_to_ic[j];
            size_t c_a          = grid.c_get_flatten_index(c_j, c_k);

            auto fnc_external = [&](double xi) -> double {
               auto fnc_internal_1 = [&](double eta) -> double {
                  auto rhophi = RnC::from_x123_to_rhophi(eta, xi - eta, -xi);
                  return -grid._dw_dx3[index](rhophi) / eta;
               };

               return integrator::integrate(fnc_internal_1, xBj, xi, int_e_r, int_e_a) / xi;
            };
            const double tmp = integrator::integrate(fnc_external, xi_min, xi_max, int_e_r, int_e_a);
            // Works both if key does and does not exist.
            std::unique_lock<std::mutex> lock(done_mutex);
            weights[c_a] += tmp;
            lock.unlock();
         });
      }
   }

   // Sector 3, only third integral (the one with switched arguments)
   for (size_t k = grid.grid_angle._delim_indexes[3]; k < grid.grid_angle._delim_indexes[4]; k++) {
      for (size_t j = 0; j < grid.grid_radius.size; j++) {

         th.AddTask([&, j, k](void) -> void {
            auto r_supp         = grid.grid_radius.get_support_weight_aj(j);
            const double xi_min = std::max(grid.grid_radius._d_info.to_phys_space(r_supp.first), xBj);
            const double xi_max = std::min(grid.grid_radius._d_info.to_phys_space(r_supp.second), 1.0);
            auto index          = grid.get_flatten_index(j, k);
            size_t c_k          = grid.grid_angle._from_iw_to_ic[k];
            size_t c_j          = grid.grid_radius._from_iw_to_ic[j];
            size_t c_a          = grid.c_get_flatten_index(c_j, c_k);

            auto fnc_external = [&](double xi) -> double {
               auto fnc_internal_1 = [&](double eta) -> double {
                  // Switched arguments, and mind the global `-` sign
                  auto rhophi = RnC::from_x123_to_rhophi(eta, -eta - xi, xi);
                  return -(grid._dw_dx3[index](rhophi)) / (eta + xi);
               };

               return (integrator::integrate(fnc_internal_1, -xBj, 0, int_e_r, int_e_a)) / xi;
            };
            const double tmp = integrator::integrate(fnc_external, xi_min, xi_max, int_e_r, int_e_a);
            // Works both if key does and does not exist.
            std::unique_lock<std::mutex> lock(done_mutex);
            weights[c_a] += tmp;
            lock.unlock();
         });
      }
   }

   // Sector 4, second and third integrals, written as integrals over t
   for (size_t k = grid.grid_angle._delim_indexes[4]; k < grid.grid_angle._delim_indexes[5]; k++) {
      for (size_t j = 0; j < grid.grid_radius.size; j++) {
         th.AddTask([&, j, k](void) -> void {
            auto r_supp         = grid.grid_radius.get_support_weight_aj(j);
            const double xi_min = std::max(grid.grid_radius._d_info.to_phys_space(r_supp.first), xBj);
            const double xi_max = std::min(grid.grid_radius._d_info.to_phys_space(r_supp.second), 1.0);
            auto index          = grid.get_flatten_index(j, k);
            size_t c_k          = grid.grid_angle._from_iw_to_ic[k];
            size_t c_j          = grid.grid_radius._from_iw_to_ic[j];
            size_t c_a          = grid.c_get_flatten_index(c_j, c_k);

            auto fnc_external = [&](double xi) -> double {
               auto fnc_internal_1 = [&](double eta) -> double {
                  auto rhophi_1 = RnC::from_x123_to_rhophi(eta, -xi, xi - eta);
                  auto rhophi_2 = RnC::from_x123_to_rhophi(xi - eta, -xi, eta);
                  return (grid._dw_dx3[index](rhophi_1) - grid._dw_dx3[index](rhophi_2)) / eta;
               };

               // yes, |t| in denominator
               return (integrator::integrate(fnc_internal_1, xBj, xi, int_e_r, int_e_a)) / xi;
            };
            const double tmp = integrator::integrate(fnc_external, xi_min, xi_max, int_e_r, int_e_a);
            // Works both if key does and does not exist.
            std::unique_lock<std::mutex> lock(done_mutex);
            weights[c_a] += tmp;
            lock.unlock();
         });
      }
   }

   // Sector 5, both first and second integrals
   for (size_t k = grid.grid_angle._delim_indexes[5]; k < grid.grid_angle._delim_indexes[6]; k++) {
      for (size_t j = 0; j < grid.grid_radius.size; j++) {
         th.AddTask([&, j, k](void) -> void {
            auto r_supp         = grid.grid_radius.get_support_weight_aj(j);
            const double xi_min = std::max(grid.grid_radius._d_info.to_phys_space(r_supp.first), xBj);
            const double xi_max = std::min(grid.grid_radius._d_info.to_phys_space(r_supp.second), 1.0);
            auto index          = grid.get_flatten_index(j, k);
            size_t c_k          = grid.grid_angle._from_iw_to_ic[k];
            size_t c_j          = grid.grid_radius._from_iw_to_ic[j];
            size_t c_a          = grid.c_get_flatten_index(c_j, c_k);

            auto fnc_external = [&](double xi) -> double {
               auto fnc_internal_1 = [&](double eta) -> double {
                  auto rhophi = RnC::from_x123_to_rhophi(xi, -eta - xi, eta);
                  return grid._dw_dx3[index](rhophi) / eta;
               };
               auto fnc_internal_2 = [&](double eta) -> double {
                  auto rhophi = RnC::from_x123_to_rhophi(xi, -eta - xi, eta);
                  return (grid._dw_dx3[index](rhophi)) / (eta + xi);
               };

               // Combines pieces from first and second integral
               return (integrator::integrate(fnc_internal_1, -xi, -xBj, int_e_r, int_e_a)
                       + integrator::integrate(fnc_internal_2, -xBj, 0, int_e_r, int_e_a))
                    / xi;
            };
            const double tmp = integrator::integrate(fnc_external, xi_min, xi_max, int_e_r, int_e_a);
            // Works both if key does and does not exist.
            std::unique_lock<std::mutex> lock(done_mutex);
            weights[c_a] += tmp;
            lock.unlock();
         });
      }
   }

   th.WaitOnJobs();

   // ! Extra factor 1/2 from convention
   for (long int i = 0; i < grid.c_size_li; i++)
      weights[i] *= 0.5;

   return weights;
}

G2Weights::G2Weights(double _xBj, const Grid2D &_grid, bool to_compute, double int_e_r, double int_e_a)
    : xBj(_xBj), grid(_grid)
{
   if (to_compute) weights = G2Weights::GetWeights(_xBj, _grid, int_e_r, int_e_a);
}

//================================================================================

D2Weights::D2Weights(const Grid2D &_grid, bool to_compute, double int_e_r, double int_e_a) : grid(_grid)
{

   if (to_compute) GetWeights(int_e_r, int_e_a);
}

void D2Weights::GetWeights(double int_e_r, double int_e_a)
{
   if (!grid.is_compliant) {
      logger(Logger::ERROR, "D2Weights: The Grid2D is not compliant to the specs. Please use only grids "
                            "generated via `generate_compliant_Grid2D`.");
   }

   const double r0 = grid.grid_radius._d_info.intervals_phys[0].first;

   using integrator            = Honeycomb::Integration::GaussKronrod<Honeycomb::Integration::GK_61>;
   unsigned int num_av_threads = std::thread::hardware_concurrency();
   unsigned int num_threads    = num_av_threads <= 2 ? 1 : num_av_threads - 2;
   ThreadPool th(num_threads);

   std::mutex done_mutex;

   weights = Eigen::VectorXd::Zero(grid.c_size_li);

   const std::vector<std::function<double(const RnC::Pair &rhophi)>> &_w_fnc = grid._w_extrap;

   //  Sector 0,
   for (size_t k = grid.grid_angle._delim_indexes[0]; k < grid.grid_angle._delim_indexes[1]; k++) {
      for (size_t j = 0; j < grid.grid_radius.size; j++) {
         th.AddTask([&, j, k](void) -> void {
            auto r_supp = grid.grid_radius.get_support_weight_aj(j);

            double lower = grid.grid_radius._d_info.to_phys_space(r_supp.first);
            if (is_near(lower, r0)) lower = 0.0;

            const double x3_r_min = std::max(lower, 0.0);
            const double x3_r_max = std::min(grid.grid_radius._d_info.to_phys_space(r_supp.second), 1.0);
            auto index            = grid.get_flatten_index(j, k);
            size_t c_k            = grid.grid_angle._from_iw_to_ic[k];
            size_t c_j            = grid.grid_radius._from_iw_to_ic[j];
            size_t c_a            = grid.c_get_flatten_index(c_j, c_k);

            auto fnc_external = [&](double x3_r) -> double {
               auto fnc_internal_1 = [&](double x1) -> double {
                  auto rhophi = RnC::from_x123_to_rhophi(x1, x3_r - x1, -x3_r);
                  return 0.25 * _w_fnc[index](rhophi);
               };

               return integrator::integrate(fnc_internal_1, 0.0, x3_r, int_e_r, int_e_a);
            };
            const double tmp = integrator::integrate(fnc_external, x3_r_min, x3_r_max, int_e_r, int_e_a);

            // Works both if key does and does not exist.
            std::unique_lock<std::mutex> lock(done_mutex);
            weights[c_a] += tmp;
            lock.unlock();
         });
      }
   }

   //  Sector 3,
   for (size_t k = grid.grid_angle._delim_indexes[3]; k < grid.grid_angle._delim_indexes[4]; k++) {
      for (size_t j = 0; j < grid.grid_radius.size; j++) {
         th.AddTask([&, j, k](void) -> void {
            auto r_supp = grid.grid_radius.get_support_weight_aj(j);

            double lower = grid.grid_radius._d_info.to_phys_space(r_supp.first);
            if (is_near(lower, r0)) lower = 0.0;

            const double x3_min = std::max(lower, 0.0);
            const double x3_max = std::min(grid.grid_radius._d_info.to_phys_space(r_supp.second), 1.0);
            auto index          = grid.get_flatten_index(j, k);
            size_t c_k          = grid.grid_angle._from_iw_to_ic[k];
            size_t c_j          = grid.grid_radius._from_iw_to_ic[j];
            size_t c_a          = grid.c_get_flatten_index(c_j, c_k);

            auto fnc_external = [&](double x3) -> double {
               auto fnc_internal_1 = [&](double x1) -> double {
                  auto rhophi = RnC::from_x123_to_rhophi(x1, -x3 - x1, x3);
                  return 0.25 * _w_fnc[index](rhophi);
               };

               return integrator::integrate(fnc_internal_1, -x3, 0.0, int_e_r, int_e_a);
            };
            const double tmp = integrator::integrate(fnc_external, x3_min, x3_max, int_e_r, int_e_a);
            // Works both if key does and does not exist.
            std::unique_lock<std::mutex> lock(done_mutex);
            weights[c_a] += tmp;
            lock.unlock();
         });
      }
   }

   //  Sector 1,
   for (size_t k = grid.grid_angle._delim_indexes[1]; k < grid.grid_angle._delim_indexes[2]; k++) {
      for (size_t j = 0; j < grid.grid_radius.size; j++) {
         th.AddTask([&, j, k](void) -> void {
            auto r_supp = grid.grid_radius.get_support_weight_aj(j);

            double lower = grid.grid_radius._d_info.to_phys_space(r_supp.first);
            if (is_near(lower, r0)) lower = 0.0;

            const double x2_min = std::max(lower, 0.0);
            const double x2_max = std::min(grid.grid_radius._d_info.to_phys_space(r_supp.second), 1.0);
            auto index          = grid.get_flatten_index(j, k);
            size_t c_k          = grid.grid_angle._from_iw_to_ic[k];
            size_t c_j          = grid.grid_radius._from_iw_to_ic[j];
            size_t c_a          = grid.c_get_flatten_index(c_j, c_k);

            auto fnc_external = [&](double x2) -> double {
               auto fnc_internal_1 = [&](double x1) -> double {
                  auto rhophi = RnC::from_x123_to_rhophi(x1, x2, -x1 - x2);
                  return 0.25 * _w_fnc[index](rhophi);
               };

               return integrator::integrate(fnc_internal_1, -x2, 0.0, int_e_r, int_e_a);
            };
            const double tmp = integrator::integrate(fnc_external, x2_min, x2_max, int_e_r, int_e_a);
            // Works both if key does and does not exist.
            std::unique_lock<std::mutex> lock(done_mutex);
            weights[c_a] += tmp;
            lock.unlock();
         });
      }
   }

   //  Sector 4,
   for (size_t k = grid.grid_angle._delim_indexes[4]; k < grid.grid_angle._delim_indexes[5]; k++) {
      for (size_t j = 0; j < grid.grid_radius.size; j++) {
         th.AddTask([&, j, k](void) -> void {
            auto r_supp = grid.grid_radius.get_support_weight_aj(j);

            double lower = grid.grid_radius._d_info.to_phys_space(r_supp.first);
            if (is_near(lower, r0)) lower = 0.0;

            const double x2_r_min = std::max(lower, 0.0);
            const double x2_r_max = std::min(grid.grid_radius._d_info.to_phys_space(r_supp.second), 1.0);
            auto index            = grid.get_flatten_index(j, k);
            size_t c_k            = grid.grid_angle._from_iw_to_ic[k];
            size_t c_j            = grid.grid_radius._from_iw_to_ic[j];
            size_t c_a            = grid.c_get_flatten_index(c_j, c_k);

            auto fnc_external = [&](double x2_r) -> double {
               auto fnc_internal_1 = [&](double x1) -> double {
                  auto rhophi = RnC::from_x123_to_rhophi(x1, -x2_r, -x1 + x2_r);
                  return 0.25 * _w_fnc[index](rhophi);
               };

               return integrator::integrate(fnc_internal_1, 0.0, x2_r, int_e_r, int_e_a);
            };
            const double tmp = integrator::integrate(fnc_external, x2_r_min, x2_r_max, int_e_r, int_e_a);
            // Works both if key does and does not exist.
            std::unique_lock<std::mutex> lock(done_mutex);
            weights[c_a] += tmp;
            lock.unlock();
         });
      }
   }

   //  Sector 2,
   for (size_t k = grid.grid_angle._delim_indexes[2]; k < grid.grid_angle._delim_indexes[3]; k++) {
      for (size_t j = 0; j < grid.grid_radius.size; j++) {
         th.AddTask([&, j, k](void) -> void {
            auto r_supp = grid.grid_radius.get_support_weight_aj(j);

            double lower = grid.grid_radius._d_info.to_phys_space(r_supp.first);
            if (is_near(lower, r0)) lower = 0.0;

            const double x1_r_min = std::max(lower, 0.0);
            const double x1_r_max = std::min(grid.grid_radius._d_info.to_phys_space(r_supp.second), 1.0);
            auto index            = grid.get_flatten_index(j, k);
            size_t c_k            = grid.grid_angle._from_iw_to_ic[k];
            size_t c_j            = grid.grid_radius._from_iw_to_ic[j];
            size_t c_a            = grid.c_get_flatten_index(c_j, c_k);

            auto fnc_external = [&](double x1_r) -> double {
               auto fnc_internal_1 = [&](double x2) -> double {
                  auto rhophi = RnC::from_x123_to_rhophi(-x1_r, x2, x1_r - x2);
                  return 0.25 * _w_fnc[index](rhophi);
               };

               return integrator::integrate(fnc_internal_1, 0.0, x1_r, int_e_r, int_e_a);
            };
            const double tmp = integrator::integrate(fnc_external, x1_r_min, x1_r_max, int_e_r, int_e_a);
            // Works both if key does and does not exist.
            std::unique_lock<std::mutex> lock(done_mutex);
            weights[c_a] += tmp;
            lock.unlock();
         });
      }
   }

   //  Sector 5,
   for (size_t k = grid.grid_angle._delim_indexes[5]; k < grid.grid_angle._delim_indexes[6]; k++) {
      for (size_t j = 0; j < grid.grid_radius.size; j++) {
         th.AddTask([&, j, k](void) -> void {
            auto r_supp = grid.grid_radius.get_support_weight_aj(j);

            double lower = grid.grid_radius._d_info.to_phys_space(r_supp.first);
            if (is_near(lower, r0)) lower = 0.0;

            const double x1_min = std::max(lower, 0.0);
            const double x1_max = std::min(grid.grid_radius._d_info.to_phys_space(r_supp.second), 1.0);
            auto index          = grid.get_flatten_index(j, k);
            size_t c_k          = grid.grid_angle._from_iw_to_ic[k];
            size_t c_j          = grid.grid_radius._from_iw_to_ic[j];
            size_t c_a          = grid.c_get_flatten_index(c_j, c_k);

            auto fnc_external = [&](double x1) -> double {
               auto fnc_internal_1 = [&](double x2) -> double {
                  auto rhophi = RnC::from_x123_to_rhophi(x1, x2, -x1 - x2);
                  return 0.25 * _w_fnc[index](rhophi);
               };

               return integrator::integrate(fnc_internal_1, -x1, 0.0, int_e_r, int_e_a);
            };
            const double tmp = integrator::integrate(fnc_external, x1_min, x1_max, int_e_r, int_e_a);
            // Works both if key does and does not exist.
            std::unique_lock<std::mutex> lock(done_mutex);
            weights[c_a] += tmp;
            lock.unlock();
         });
      }
   }

   th.WaitOnJobs();
}

double D2Weights::ComputeSingleQuark(const Eigen::VectorXd &_f) const
{
   return weights.dot(_f);
}

//================================================================================

D2WeightsCutted::D2WeightsCutted(const Grid2D &_grid, bool to_compute, double int_e_r, double int_e_a)
    : grid(_grid)
{

   center_approx = 0;
   if (to_compute) GetWeights(int_e_r, int_e_a);
}

void D2WeightsCutted::GetWeights(double int_e_r, double int_e_a)
{
   if (!grid.is_compliant) {
      logger(Logger::ERROR,
             "D2WeightsCutted: The Grid2D is not compliant to the specs. Please use only grids "
             "generated via `generate_compliant_Grid2D`.");
   }
   const double r0 = grid.grid_radius._d_info.intervals_phys[0].first;

   center_approx = 3.0 * r0 * r0 / 8.0;

   using integrator            = Honeycomb::Integration::GaussKronrod<Honeycomb::Integration::GK_61>;
   unsigned int num_av_threads = std::thread::hardware_concurrency();
   unsigned int num_threads    = num_av_threads <= 2 ? 1 : num_av_threads - 2;
   ThreadPool th(num_threads);

   std::mutex done_mutex;

   weights = Eigen::VectorXd::Zero(grid.c_size_li);

   const std::vector<std::function<double(const RnC::Pair &rhophi)>> &_w_fnc = grid._w;

   //  Sector 0,
   for (size_t k = grid.grid_angle._delim_indexes[0]; k < grid.grid_angle._delim_indexes[1]; k++) {
      for (size_t j = 0; j < grid.grid_radius.size; j++) {
         th.AddTask([&, j, k](void) -> void {
            auto r_supp = grid.grid_radius.get_support_weight_aj(j);

            double lower = grid.grid_radius._d_info.to_phys_space(r_supp.first);

            const double x3_r_min = std::max(lower, 0.0);
            const double x3_r_max = std::min(grid.grid_radius._d_info.to_phys_space(r_supp.second), 1.0);
            auto index            = grid.get_flatten_index(j, k);
            size_t c_k            = grid.grid_angle._from_iw_to_ic[k];
            size_t c_j            = grid.grid_radius._from_iw_to_ic[j];
            size_t c_a            = grid.c_get_flatten_index(c_j, c_k);

            auto fnc_external = [&](double x3_r) -> double {
               auto fnc_internal_1 = [&](double x1) -> double {
                  auto rhophi = RnC::from_x123_to_rhophi(x1, x3_r - x1, -x3_r);
                  return 0.25 * _w_fnc[index](rhophi);
               };

               return integrator::integrate(fnc_internal_1, 0.0, x3_r, int_e_r, int_e_a);
            };
            const double tmp = integrator::integrate(fnc_external, x3_r_min, x3_r_max, int_e_r, int_e_a);
            // Works both if key does and does not exist.
            std::unique_lock<std::mutex> lock(done_mutex);
            weights[c_a] += tmp;
            lock.unlock();
         });
      }
   }

   //  Sector 3,
   for (size_t k = grid.grid_angle._delim_indexes[3]; k < grid.grid_angle._delim_indexes[4]; k++) {
      for (size_t j = 0; j < grid.grid_radius.size; j++) {
         th.AddTask([&, j, k](void) -> void {
            auto r_supp = grid.grid_radius.get_support_weight_aj(j);

            double lower = grid.grid_radius._d_info.to_phys_space(r_supp.first);

            const double x3_min = std::max(lower, 0.0);
            const double x3_max = std::min(grid.grid_radius._d_info.to_phys_space(r_supp.second), 1.0);
            auto index          = grid.get_flatten_index(j, k);
            size_t c_k          = grid.grid_angle._from_iw_to_ic[k];
            size_t c_j          = grid.grid_radius._from_iw_to_ic[j];
            size_t c_a          = grid.c_get_flatten_index(c_j, c_k);

            auto fnc_external = [&](double x3) -> double {
               auto fnc_internal_1 = [&](double x1) -> double {
                  auto rhophi = RnC::from_x123_to_rhophi(x1, -x3 - x1, x3);
                  return 0.25 * _w_fnc[index](rhophi);
               };

               return integrator::integrate(fnc_internal_1, -x3, 0.0, int_e_r, int_e_a);
            };
            const double tmp = integrator::integrate(fnc_external, x3_min, x3_max, int_e_r, int_e_a);
            // Works both if key does and does not exist.
            std::unique_lock<std::mutex> lock(done_mutex);
            weights[c_a] += tmp;
            lock.unlock();
         });
      }
   }

   //  Sector 1,
   for (size_t k = grid.grid_angle._delim_indexes[1]; k < grid.grid_angle._delim_indexes[2]; k++) {
      for (size_t j = 0; j < grid.grid_radius.size; j++) {
         th.AddTask([&, j, k](void) -> void {
            auto r_supp = grid.grid_radius.get_support_weight_aj(j);

            double lower = grid.grid_radius._d_info.to_phys_space(r_supp.first);

            const double x2_min = std::max(lower, 0.0);
            const double x2_max = std::min(grid.grid_radius._d_info.to_phys_space(r_supp.second), 1.0);
            auto index          = grid.get_flatten_index(j, k);
            size_t c_k          = grid.grid_angle._from_iw_to_ic[k];
            size_t c_j          = grid.grid_radius._from_iw_to_ic[j];
            size_t c_a          = grid.c_get_flatten_index(c_j, c_k);

            auto fnc_external = [&](double x2) -> double {
               auto fnc_internal_1 = [&](double x1) -> double {
                  auto rhophi = RnC::from_x123_to_rhophi(x1, x2, -x1 - x2);
                  return 0.25 * _w_fnc[index](rhophi);
               };

               return integrator::integrate(fnc_internal_1, -x2, 0.0, int_e_r, int_e_a);
            };
            const double tmp = integrator::integrate(fnc_external, x2_min, x2_max, int_e_r, int_e_a);
            // Works both if key does and does not exist.
            std::unique_lock<std::mutex> lock(done_mutex);
            weights[c_a] += tmp;
            lock.unlock();
         });
      }
   }

   //  Sector 4,
   for (size_t k = grid.grid_angle._delim_indexes[4]; k < grid.grid_angle._delim_indexes[5]; k++) {
      for (size_t j = 0; j < grid.grid_radius.size; j++) {
         th.AddTask([&, j, k](void) -> void {
            auto r_supp = grid.grid_radius.get_support_weight_aj(j);

            double lower = grid.grid_radius._d_info.to_phys_space(r_supp.first);

            const double x2_r_min = std::max(lower, 0.0);
            const double x2_r_max = std::min(grid.grid_radius._d_info.to_phys_space(r_supp.second), 1.0);
            auto index            = grid.get_flatten_index(j, k);
            size_t c_k            = grid.grid_angle._from_iw_to_ic[k];
            size_t c_j            = grid.grid_radius._from_iw_to_ic[j];
            size_t c_a            = grid.c_get_flatten_index(c_j, c_k);

            auto fnc_external = [&](double x2_r) -> double {
               auto fnc_internal_1 = [&](double x1) -> double {
                  auto rhophi = RnC::from_x123_to_rhophi(x1, -x2_r, -x1 + x2_r);
                  return 0.25 * _w_fnc[index](rhophi);
               };

               return integrator::integrate(fnc_internal_1, 0.0, x2_r, int_e_r, int_e_a);
            };
            const double tmp = integrator::integrate(fnc_external, x2_r_min, x2_r_max, int_e_r, int_e_a);
            // Works both if key does and does not exist.
            std::unique_lock<std::mutex> lock(done_mutex);
            weights[c_a] += tmp;
            lock.unlock();
         });
      }
   }

   //  Sector 2,
   for (size_t k = grid.grid_angle._delim_indexes[2]; k < grid.grid_angle._delim_indexes[3]; k++) {
      for (size_t j = 0; j < grid.grid_radius.size; j++) {
         th.AddTask([&, j, k](void) -> void {
            auto r_supp = grid.grid_radius.get_support_weight_aj(j);

            double lower = grid.grid_radius._d_info.to_phys_space(r_supp.first);

            const double x1_r_min = std::max(lower, 0.0);
            const double x1_r_max = std::min(grid.grid_radius._d_info.to_phys_space(r_supp.second), 1.0);
            auto index            = grid.get_flatten_index(j, k);
            size_t c_k            = grid.grid_angle._from_iw_to_ic[k];
            size_t c_j            = grid.grid_radius._from_iw_to_ic[j];
            size_t c_a            = grid.c_get_flatten_index(c_j, c_k);

            auto fnc_external = [&](double x1_r) -> double {
               auto fnc_internal_1 = [&](double x2) -> double {
                  auto rhophi = RnC::from_x123_to_rhophi(-x1_r, x2, x1_r - x2);
                  return 0.25 * _w_fnc[index](rhophi);
               };

               return integrator::integrate(fnc_internal_1, 0.0, x1_r, int_e_r, int_e_a);
            };
            const double tmp = integrator::integrate(fnc_external, x1_r_min, x1_r_max, int_e_r, int_e_a);
            // Works both if key does and does not exist.
            std::unique_lock<std::mutex> lock(done_mutex);
            weights[c_a] += tmp;
            lock.unlock();
         });
      }
   }

   //  Sector 5,
   for (size_t k = grid.grid_angle._delim_indexes[5]; k < grid.grid_angle._delim_indexes[6]; k++) {
      for (size_t j = 0; j < grid.grid_radius.size; j++) {
         th.AddTask([&, j, k](void) -> void {
            auto r_supp = grid.grid_radius.get_support_weight_aj(j);

            double lower = grid.grid_radius._d_info.to_phys_space(r_supp.first);

            const double x1_min = std::max(lower, 0.0);
            const double x1_max = std::min(grid.grid_radius._d_info.to_phys_space(r_supp.second), 1.0);
            auto index          = grid.get_flatten_index(j, k);
            size_t c_k          = grid.grid_angle._from_iw_to_ic[k];
            size_t c_j          = grid.grid_radius._from_iw_to_ic[j];
            size_t c_a          = grid.c_get_flatten_index(c_j, c_k);

            auto fnc_external = [&](double x1) -> double {
               auto fnc_internal_1 = [&](double x2) -> double {
                  auto rhophi = RnC::from_x123_to_rhophi(x1, x2, -x1 - x2);
                  return 0.25 * _w_fnc[index](rhophi);
               };

               return integrator::integrate(fnc_internal_1, -x1, 0.0, int_e_r, int_e_a);
            };
            const double tmp = integrator::integrate(fnc_external, x1_min, x1_max, int_e_r, int_e_a);
            // Works both if key does and does not exist.
            std::unique_lock<std::mutex> lock(done_mutex);
            weights[c_a] += tmp;
            lock.unlock();
         });
      }
   }

   th.WaitOnJobs();
}

double D2WeightsCutted::ComputeSingleQuark(const Eigen::VectorXd &_f) const
{
   const double tmp = weights.dot(_f);
   double ave_r0    = 0;
   double n         = 0;
   for (size_t c_a = 0; c_a < grid.c_size; c_a++) {

      auto [c_jP, c_iP] = grid.c_get_double_index(c_a);
      // size_t c_jP   = grid.grid_radius._from_iw_to_ic[jP];

      if (c_jP != 0) continue;

      ave_r0 += _f[c_a];
      n++;
   }
   ave_r0 /= n;

   return ave_r0 * center_approx + tmp;
}

double D2WeightsCutted::ComputeSingleQuark_NoCorrections(const Eigen::VectorXd &_f) const
{
   return weights.dot(_f);
}

//================================================================================

ELTWeights::ELTWeights(const Grid2D &_grid, bool to_compute, double int_e_r, double int_e_a) : grid(_grid)
{

   center_approx = 0;
   if (to_compute) GetWeights(int_e_r, int_e_a);
}

void ELTWeights::GetWeights(double int_e_r, double int_e_a)
{
   if (!grid.is_compliant) {
      logger(Logger::ERROR, "ELTWeights: The Grid2D is not compliant to the specs. Please use only grids "
                            "generated via `generate_compliant_Grid2D`.");
   }
   const double r0 = grid.grid_radius._d_info.intervals_phys[0].first;
   center_approx   = r0 / 2.0;

   using integrator            = Honeycomb::Integration::GaussKronrod<Honeycomb::Integration::GK_61>;
   unsigned int num_av_threads = std::thread::hardware_concurrency();
   unsigned int num_threads    = num_av_threads <= 2 ? 1 : num_av_threads - 2;
   ThreadPool th(num_threads);

   std::mutex done_mutex;

   weights = Eigen::VectorXd::Zero(grid.c_size_li);

   const std::vector<std::function<double(const RnC::Pair &rhophi)>> &_w_fnc = grid._w;

   //  Sector 1,
   for (size_t k = grid.grid_angle._delim_indexes[1]; k < grid.grid_angle._delim_indexes[2]; k++) {
      for (size_t j = 0; j < grid.grid_radius.size; j++) {
         th.AddTask([&, j, k](void) -> void {
            auto r_supp = grid.grid_radius.get_support_weight_aj(j);

            double lower = grid.grid_radius._d_info.to_phys_space(r_supp.first);
            // if (is_near(lower, r0)) lower = 0.0;

            const double x2_min = std::max(lower, 0.0);
            const double x2_max = std::min(grid.grid_radius._d_info.to_phys_space(r_supp.second), 1.0);
            auto index          = grid.get_flatten_index(j, k);
            size_t c_k          = grid.grid_angle._from_iw_to_ic[k];
            size_t c_j          = grid.grid_radius._from_iw_to_ic[j];
            size_t c_a          = grid.c_get_flatten_index(c_j, c_k);

            auto fnc_external = [&](double x2) -> double {
               auto fnc_internal_1 = [&](double x1) -> double {
                  auto rhophi = RnC::from_x123_to_rhophi(x1, x2, -x1 - x2);
                  return (-x1 / sq(x2)) * _w_fnc[index](rhophi);
               };

               return integrator::integrate(fnc_internal_1, -x2, 0.0, int_e_r, int_e_a);
            };
            const double tmp = integrator::integrate(fnc_external, x2_min, x2_max, int_e_r, int_e_a);
            // Works both if key does and does not exist.
            std::unique_lock<std::mutex> lock(done_mutex);
            weights[c_a] += tmp;
            lock.unlock();
         });
      }
   }

   th.WaitOnJobs();
}

double ELTWeights::ComputeSingleQuark_NoCorrections(const Eigen::VectorXd &_f) const
{
   return weights.dot(_f);
}

double ELTWeights::ComputeSingleQuark(const Eigen::VectorXd &_f) const
{
   const double tmp = weights.dot(_f);
   double ave_r0    = 0;
   double n         = 0;
   for (size_t k = grid.grid_angle._delim_indexes[1]; k < grid.grid_angle._delim_indexes[2]; k++) {

      size_t c_k = grid.grid_angle._from_iw_to_ic[k];
      size_t c_a = grid.c_get_flatten_index(0, c_k);

      ave_r0 += _f[c_a];
      n++;
   }
   ave_r0 /= n;

   return ave_r0 * center_approx + tmp;
}

D2WeightsPartialIntegral::D2WeightsPartialIntegral(const Grid2D &_grid, double _a, double _b, bool to_compute,
                                                   double int_e_r, double int_e_a)
    : grid(_grid), a(_a), b(_b)
{

   if (to_compute) GetWeights(int_e_r, int_e_a);
}

void D2WeightsPartialIntegral::GetWeights(double int_e_r, double int_e_a)
{
   // \int_a^b dx 3 x^2 w_g2(x)
   const double center      = (b + a) / 2.0;
   const double half_length = (b / a) / 2.0;
   const size_t n           = Integration::GK_41::size;

   weights = Integration::GK_41::w_gk[n - 1] * 3.0 * center * center
           * G2Weights::GetWeights(center, grid, int_e_r, int_e_a);

   for (size_t i = 0; i < n - 1; i++) {
      const double x_m    = center - half_length * Integration::GK_41::x_gk[i];
      const double x_p    = center + half_length * Integration::GK_41::x_gk[i];
      const double pref_m = 3.0 * x_m * x_m * Integration::GK_41::w_gk[i];
      const double pref_p = 3.0 * x_p * x_p * Integration::GK_41::w_gk[i];

      weights += pref_m * G2Weights::GetWeights(x_m, grid, int_e_r, int_e_a);
      weights += pref_p * G2Weights::GetWeights(x_p, grid, int_e_r, int_e_a);
   }

   weights *= half_length;
}

double D2WeightsPartialIntegral::ComputeSingleQuark(const Eigen::VectorXd &_f) const
{
   return weights.dot(_f);
}

} // namespace Honeycomb