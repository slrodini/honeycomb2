#ifndef G2_WEIGHTS_HPP
#define G2_WEIGHTS_HPP

#include <honeycomb2/discretization.hpp>
#include <honeycomb2/gauss_kronrod.hpp>

namespace Honeycomb
{
struct G2Weights {

   G2Weights(double _xBj, const Grid2D &_grid, double int_e_r = 1.0e-8, double int_e_a = 1.0e-8);

   const double xBj;
   const Grid2D &grid;
   Eigen::VectorXd weights;
};

} // namespace Honeycomb

#endif
