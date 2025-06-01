#ifndef HC2_FOREIGN_INTERFACE
#define HC2_FOREIGN_INTERFACE

#include <cmath>
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
   ForeignInterfaceState()
       : discr(grid), d2_weights(grid, false), elt_weights(grid, false), _unloaded(false) {};

   void Evolve();
   double GetDistribution(OutputModel::FNC f, double Q2, double x1, double x2, double x3);
   void Unload();

   Grid2D grid;
   Discretization discr;
   InputModel model;
   std::function<double(double)> as_fnc = [](double) -> double {
      logger(Logger::ERROR, "Alpha_s has not been correctly provided!");
      return NAN;
   };

   std::vector<double> thresholds = {0, 0, 0, 1.6129, 17.4724, 1.0e+6};
   D2Weights d2_weights;
   ELTWeights elt_weights;
   std::vector<double> interm_scales; // scales, in GeV^2 (i.e. Q0^2, Qf^2)
   // These are for staggered evolution:
   // Q0^2 -> Q1^2 -> Q2^2 -> .. -> Qf^2
   std::vector<EvOp> evol_op;
   std::vector<Solution> _solutions; // For any number of scales.

   // Computed
   size_t nf_initial_scale;

   // Only one model is accepted for the moment.
   InputModel _models;
   std::vector<OutputModel> _fin_models;

   bool _unloaded = false;
};
} // namespace Honeycomb

#endif