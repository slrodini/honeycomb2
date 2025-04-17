#include <honeycomb2/g2_weights.hpp>
#include <honeycomb2/thread_pool.hpp>
#include <honeycomb2/utilities.hpp>

namespace Honeycomb
{

G2Weights::G2Weights(double _xBj, const Grid2D &_grid, double int_e_r, double int_e_a)
    : xBj(_xBj), grid(_grid)
{

   if (!_grid.is_compliant) {
      logger(Logger::ERROR, "G2Weights: The Grid2D is not compliant to the specs. Please use only grids "
                            "generated via `generate_compliant_Grid2D`.");
   }
   using integrator = Honeycomb::Integration::GaussKronrod<Honeycomb::Integration::GK_61>;

   ThreadPool th(10);

   std::mutex done_mutex;

   weights = Eigen::VectorXd::Zero(grid.c_size_li);

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
}
} // namespace Honeycomb