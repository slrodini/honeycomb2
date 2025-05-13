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

} // namespace Honeycomb

#endif
