#pragma once

#include <honeycomb2/discretization.hpp>
#include <honeycomb2/kernel_functions.hpp>
#include <honeycomb2/cereal_extension.hpp>

#include <Eigen/Core>

namespace Honeycomb
{

template <typename T>
concept Arithmetic = std::is_arithmetic_v<T>;

// NOTE: Temporary struct, for laying down the structure of the program.
// TODO: Update to correct solution
struct Solution {

   const Discretization *discretization;
   size_t nf;
   std::map<int, Eigen::VectorXd> _distr;

   Solution(const Discretization *_discr, std::map<int, std::function<double(double, double, double)>> models,
            size_t nf)
       : discretization(_discr), nf(nf)
   {
      _distr[0] = _discr->discretize(models.at(0));
   };

   void PushNewFlavor()
   {
   }

   // Right-hand multiplication (Solution * scalar)
   template <Arithmetic T>
   Solution operator*(T scalar) const
   {
      Solution result(*this);
      // result._distr._fj = scalar * _distr._fj;
      return result;
   }

   // Right-hand division (Solution / scalar)
   template <Arithmetic T>
   Solution operator/(T scalar) const
   {
      Solution result(*this);
      // result._distr._fj = _distr._fj / scalar;
      return result;
   }

   // Left-hand multiplication (scalar * Solution)
   template <Arithmetic T>
   friend Solution operator*(T scalar, const Solution &obj)
   {
      Solution result(obj);
      // result._distr._fj = scalar * obj._distr._fj;
      return result;
   }

   // First, define the addition operator
   Solution operator+(const Solution &other) const
   {
      Solution result(*this);
      // result._distr._fj = other._distr._fj + _distr._fj;
      return result;
   }

   // Then, define the compound assignment operator
   Solution &operator+=(const Solution &other)
   {
      // _distr._fj += other._distr._fj;
      return *this;
   }
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