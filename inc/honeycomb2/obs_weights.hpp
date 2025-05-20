#ifndef G2_WEIGHTS_HPP
#define G2_WEIGHTS_HPP

#include <honeycomb2/discretization.hpp>
#include <honeycomb2/solution.hpp>
#include <honeycomb2/gauss_kronrod.hpp>

namespace Honeycomb
{
struct G2Weights {

   G2Weights(double _xBj, const Grid2D &_grid, double int_e_r = 1.0e-8, double int_e_a = 1.0e-8);

   const double xBj;
   const Grid2D &grid;
   Eigen::VectorXd weights;
};

// Weights for the d_2 observable, defined as
// 3\int_0^1 dx x^2 g_2(x, Q) = \int Dx \mathfrak{S}^+
struct D2Weights {

   D2Weights(const Grid2D &_grid, double int_e_r = 1.0e-8, double int_e_a = 1.0e-8);

   double ComputeSingleQuark(const Eigen::VectorXd &_f) const;

   const Grid2D &grid;
   Eigen::VectorXd weights;
};

struct D2WeightsCutted {

   D2WeightsCutted(const Grid2D &_grid, double int_e_r = 1.0e-8, double int_e_a = 1.0e-8);

   double ComputeSingleQuark(const Eigen::VectorXd &_f) const;
   double ComputeSingleQuark_NoCorrections(const Eigen::VectorXd &_f) const;

   const Grid2D &grid;
   Eigen::VectorXd weights;
   double center_approx;
};

// Weights for  Efremov-LeaderTeryaev sum rule, defined as
// 2\int_0^1 dx x g_2(x, Q) = \int_0^1 dx_2 \int_{-x_2}^0 dx_1 (-x_1/x_2^2)\mathfrak{S}^+
struct EFTWeights {

   EFTWeights(const Grid2D &_grid, double int_e_r = 1.0e-8, double int_e_a = 1.0e-8);

   double ComputeSingleQuark(const Eigen::VectorXd &_f) const;
   double ComputeSingleQuark_NoCorrections(const Eigen::VectorXd &_f) const;

   const Grid2D &grid;
   Eigen::VectorXd weights;
   double center_approx;
};

} // namespace Honeycomb

#endif
