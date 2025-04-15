#ifndef KERNELS_HPP
#define KERNELS_HPP

#include <honeycomb2/discretization.hpp>
#include <honeycomb2/kernel_functions.hpp>
#include <honeycomb2/cereal_extension.hpp>

#include <Eigen/Core>

namespace Honeycomb
{

template <typename T>
concept Arithmetic = std::is_arithmetic_v<T>;

inline double zero_function(double a, double b, double c)
{
   return 0;
}

struct InputModel {

   enum FNC { T_UP, T_DN, T_ST, T_CH, T_BM, T_TP, DT_UP, DT_DN, DT_ST, DT_CH, DT_BM, DT_TP, T_P_GL, T_M_GL };
   void SetModel(FNC f, std::function<double(double, double, double)> model);

   std::function<double(double, double, double)> T_up = zero_function;
   std::function<double(double, double, double)> T_dn = zero_function;
   std::function<double(double, double, double)> T_st = zero_function;
   std::function<double(double, double, double)> T_ch = zero_function;
   std::function<double(double, double, double)> T_bm = zero_function;
   std::function<double(double, double, double)> T_tp = zero_function;

   std::function<double(double, double, double)> DT_up = zero_function;
   std::function<double(double, double, double)> DT_dn = zero_function;
   std::function<double(double, double, double)> DT_st = zero_function;
   std::function<double(double, double, double)> DT_ch = zero_function;
   std::function<double(double, double, double)> DT_bm = zero_function;
   std::function<double(double, double, double)> DT_tp = zero_function;

   std::function<double(double, double, double)> T_p_gl = zero_function;
   std::function<double(double, double, double)> T_m_gl = zero_function;
};

struct Kernels {

   Kernels(const Grid2D &g, double _Nc, bool to_compute = true);

   const Grid2D &grid;

   const double Nc;
   const double CA;
   const double CF;

   Eigen::MatrixXd H_NS, H_d13, H_gg_p, H_gg_m, H_qg_p, H_qg_m, H_gq_p, H_gq_m;
   Eigen::MatrixXd H_test_1, H_test_2;

   template <class Archive>
   void save(Archive &archive) const
   {
      // This saves the dimensions AND all the _x123 points, for future checks
      archive(grid);

      // Save the number of colors
      archive(Nc);

      // Explicitly archive matrix elements
      for (int i = 0; i < H_NS.rows(); ++i) {
         for (int j = 0; j < H_NS.cols(); ++j) {
            archive(H_NS(i, j));
            archive(H_d13(i, j));
            archive(H_gg_p(i, j));
            archive(H_gg_m(i, j));
            archive(H_qg_p(i, j));
            archive(H_qg_m(i, j));
            archive(H_gq_p(i, j));
            archive(H_gq_m(i, j));
         }
      }
   }

   template <class Archive>
   void load(Archive &archive)
   {
      Grid2D tmp_grid;
      // This laods the dimensions and the points
      archive(tmp_grid);
      if (tmp_grid.c_size != grid.c_size || tmp_grid.size != grid.size) {
         logger(Logger::ERROR,
                "Trying to load incompatible kernels with the given grid. The dimensions do not match. Aborting.");
      }
      for (size_t i = 0; i < grid.c_size; i++) {
         if (grid._x123[i].distance(tmp_grid._x123[i]) > 1.0e-12) {
            logger(Logger::ERROR,
                   "Trying to load incompatible kernels with the given grid. The grid points do not match. Aborting.");
         }
      }

      for (size_t i = 0; i < grid.size; i++) {
         if (grid._x123_minmax[i].first.distance(tmp_grid._x123_minmax[i].first) > 1.0e-12) {
            logger(Logger::ERROR, "Trying to load incompatible kernels with the given grid. The xmins for the weights "
                                  "do not match. Aborting.");
         }
         if (grid._x123_minmax[i].second.distance(tmp_grid._x123_minmax[i].second) > 1.0e-12) {
            logger(Logger::ERROR, "Trying to load incompatible kernels with the given grid. The xmaxs for the weights "
                                  "do not match. Aborting.");
         }
      }

      double tmp_Nc = 0;
      // Load the number of colors
      archive(tmp_Nc);

      if (std::fabs(tmp_Nc - Nc) > 1.0e-12) {
         logger(Logger::ERROR,
                "Trying to load incompatible kernels with the given grid. Numbers of color do not match. Aborting.");
      }

      // Explicitly archive matrix elements
      for (int i = 0; i < H_NS.rows(); ++i) {
         for (int j = 0; j < H_NS.cols(); ++j) {
            archive(H_NS(i, j));
            archive(H_d13(i, j));
            archive(H_gg_p(i, j));
            archive(H_gg_m(i, j));
            archive(H_qg_p(i, j));
            archive(H_qg_m(i, j));
            archive(H_gq_p(i, j));
            archive(H_gq_m(i, j));
         }
      }
   }
};

struct Solution {

   const Discretization *_discretization;
   size_t nf;

   // Gluon, Singlet, NS_1, NS_2, ...
   std::vector<Eigen::VectorXd> _distr_p;
   std::vector<Eigen::VectorXd> _distr_m;

   Solution(const Discretization *discretization, std::vector<std::function<double(double, double, double)>> models,
            size_t nf)
       : _discretization(discretization), nf(nf) {
            // _distr[0] = _discr->discretize(models.at(0));
         };

   void PushNewFlavor()
   {
   }

   // Runge-Kutta methods
   void _copy(const Solution &other)
   {
      // TODO: copy solutions
   }
   void _plus_eq(double x, const Solution &other)
   {
      // TODO: implement S += x*other
   }
   void _ker_mul(double pref, const Kernels &ker)
   {
      // TODO: implement S = pref * (ker * S);
   }
};

// The Chiral-Odd case is special, no need to keep track of anything but the kernel matrix itself
// So, I have specialized function to compute it
Eigen::MatrixXd get_CO_kernel(const Grid2D &g, double _Nc);

template <class Archive>
void save_kernels(const Kernels &k, const std::string &file_name)
{
   SaveChecksumArchive<Kernels, Archive>(k, file_name);
}

template <class Archive>
Kernels load_kernels(const std::string &file_name, const Grid2D &g, double _Nc)
{
   Kernels ker(g, _Nc, false);
   if (!LoadAndVerify<Kernels, Archive>(file_name, ker)) {
      logger(Logger::ERROR, "I was not able to correctly load the cereal archive.");
   }
   return ker;
}

} // namespace Honeycomb

#endif // KERNELS_HPP