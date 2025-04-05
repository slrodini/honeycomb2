#include <honeycomb2/kernels.hpp>
// #include <omp.h>

// TODO: The logic seems sound becuase Hplus12 convolved with a test function is prefectly smooth
// TODO: The problem is in Hhat12, which is likely the plus prescription that is giving issues.

namespace Honeycomb
{

Kernels::Kernels(const Grid2D &g, double _Nc) : grid(g), Nc(_Nc), CA(_Nc), CF((_Nc * _Nc - 1) / (2.0 * _Nc))
{

   H_NS   = Eigen::MatrixXd::Zero(g.c_size, g.c_size);
   H_CO   = Eigen::MatrixXd::Zero(g.c_size, g.c_size);
   H_d13  = Eigen::MatrixXd::Zero(g.c_size, g.c_size);
   H_gg_p = Eigen::MatrixXd::Zero(g.c_size, g.c_size);
   H_gg_m = Eigen::MatrixXd::Zero(g.c_size, g.c_size);
   H_qg_p = Eigen::MatrixXd::Zero(g.c_size, g.c_size);
   H_qg_m = Eigen::MatrixXd::Zero(g.c_size, g.c_size);
   H_gq_p = Eigen::MatrixXd::Zero(g.c_size, g.c_size);
   H_gq_m = Eigen::MatrixXd::Zero(g.c_size, g.c_size);

   // H_{aa'} = H_{[c_i c_j],[c_i' c_j']} where c_i over angles and c_j over radii
   // Then c_. can be obtained from weight index via mapping,
   // So, I loop over weight index to have simple expressions for integrals

#pragma omp parallel for
   for (size_t c_a = 0; c_a < g.c_size; c_a++) {

      for (size_t aP = 0; aP < g.size; aP++) {

         // std::fprintf(stderr, "\r\033[2K"); // Clean line
         // std::fprintf(stderr, "Doing: %zu, %zu", a, aP);

         auto [jP, iP] = g.get_double_index(aP);
         size_t c_iP   = g.grid_angle._from_iw_to_ic[iP];
         size_t c_jP   = g.grid_radius._from_iw_to_ic[jP];

         size_t c_aP = g.c_get_flatten_index(c_jP, c_iP);
         double res  = 0;

         res = Hhat12::subtracted_integrate(c_a, aP, g);
         // res = Hplus12::integrate(c_a, aP, g);
#pragma omp critical
         {
            H_NS(c_a, c_aP) += res;
         }
      }
   }
   // std::fprintf(stderr, "\n");

   // Here compute Integrate H(x123_[c_i c_j], v) \times w_{ij}(\rho, \phi)
   // transformation of w_{ij} is done in constructor.
};

} // namespace Honeycomb