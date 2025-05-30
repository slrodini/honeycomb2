#ifndef HC2_FOREIGN_INTERFACE
#define HC2_FOREIGN_INTERFACE

#include <honeycomb2/default.hpp>
#include <honeycomb2/utilities.hpp>
#include <honeycomb2/random_engine.hpp>
#include <honeycomb2/timer.hpp>
#include <honeycomb2/checksum.hpp>
#include <honeycomb2/cereal_extension.hpp>
#include <honeycomb2/discretization.hpp>

#include <honeycomb2/runge_kutta.hpp>
#include <honeycomb2/gauss_kronrod.hpp>

#include <honeycomb2/thread_pool.hpp>

#include <honeycomb2/kernels.hpp>
#include <honeycomb2/solution.hpp>

#include <honeycomb2/config_parser.hpp>

#include <honeycomb2/obs_weights.hpp>

namespace Honeycomb
{
struct ForeignInterfaceState {
   ForeignInterfaceState() : discr(grid), d2_weights(grid, false), elt_weights(grid, false) {};

   Grid2D grid;
   Discretization discr;
   InputModel model;
   std::function<double(double)> as_fnc = [](double t) -> double {
      return 1.0 / (11.0 * (t + 3.0));
   };

   std::vector<double> thresholds = {0, 0, 0, 1.6129, 17.4724, 1.0e+6};
   D2Weights d2_weights;
   ELTWeights elt_weights;
   double Q02;
   std::vector<double> fin_scales;   // final scales, in GeV^2 (i.e. Qf^2)
   std::vector<EvOp> evol_op;        // For any number of scales.
   std::vector<Solution> _solutions; // For any number of scales.
};
} // namespace Honeycomb

#endif