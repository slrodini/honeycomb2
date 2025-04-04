#pragma once

#include <honeycomb2/discretization.hpp>
#include <honeycomb2/kernel_functions.hpp>

#include <eigen3/Eigen/Dense>

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

   Kernels(const Grid2D &g, double _Nc);

   const Grid2D &grid;

   const double Nc;
   const double CA;
   const double CF;

   Eigen::MatrixXd H_NS, H_CO, H_d13, H_gg_p, H_gg_m, H_qg_p, H_qg_m, H_gq_p, H_gq_m;
};

} // namespace Honeycomb