#ifndef HC2_KERNELS_HPP
#define HC2_KERNELS_HPP

#include "honeycomb2/cereal/archives/portable_binary.hpp"
#include <honeycomb2/discretization.hpp>
#include <honeycomb2/kernel_functions.hpp>
#include <honeycomb2/cereal_extension.hpp>

#include <honeycomb2/Eigen/Core>

namespace Honeycomb
{

template <typename T>
concept Arithmetic = std::is_arithmetic_v<T>;

struct Kernels {

   Kernels(const Grid2D &g, double _Nc, bool to_compute = true);

   void ComputeKernels();
   void ComputeTestKernel();

   const Grid2D &grid;

   const double Nc;
   const double CA;
   const double CF;

   Eigen::MatrixXd TestMatrix;
   // Note: H_qg and H_d13 do not have nf, must be multiplied outside
   Eigen::MatrixXd H_NS, H_d13, H_gg_p, H_gg_m, H_qg_p, H_qg_m, H_gq_p, H_gq_m;

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

      check_grid_compatibility(tmp_grid, grid);

      double tmp_Nc = 0;
      // Load the number of colors
      archive(tmp_Nc);

      if (std::fabs(tmp_Nc - Nc) > 1.0e-12) {
         logger(Logger::ERROR, "Trying to load incompatible kernels with the given grid. Numbers of color do "
                               "not match. Aborting.");
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

struct MergedKernelsFixedNf {
   MergedKernelsFixedNf(const Kernels &ker, size_t nf);

   size_t nf;
   double beta0;
   Eigen::MatrixXd H_NS;
   Eigen::MatrixXd H_S_P, H_S_M;
};

// The Chiral-Odd case is special, no need to keep track of anything but the kernel matrix itself
// So, I have specialized function to compute it
Eigen::MatrixXd get_CO_kernel(const Grid2D &g, double _Nc);

inline void save_kernels(const Kernels &k, const std::string &file_name)
{
   SaveChecksumArchive<Kernels, cereal::PortableBinaryOutputArchive>(k, file_name);
}

inline Kernels load_kernels(const std::string &file_name, const Grid2D &g, double _Nc)
{
   Kernels ker(g, _Nc, false);
   if (!LoadAndVerify<Kernels, cereal::PortableBinaryInputArchive>(file_name, ker)) {
      logger(Logger::WARNING, "I was not able to correctly load the cereal archive " + file_name
                                  + " containing the kernels. "
                                    "Compute (and possibly overwrite)");
      ker.ComputeKernels();
      save_kernels(ker, file_name);
   }
   return ker;
}

} // namespace Honeycomb

#endif // HC2_KERNELS_HPP