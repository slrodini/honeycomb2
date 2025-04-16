#ifndef SOLUTION_HPP
#define SOLUTION_HPP

#include <honeycomb2/default.hpp>
#include <honeycomb2/kernels.hpp>
#include <honeycomb2/discretization.hpp>
#include <Eigen/Core>

namespace Honeycomb
{

inline double zero_function(double a, double b, double c)
{
   return 0;
}

struct InputModel {

   enum FNC { T_UP, T_DN, T_ST, T_CH, T_BM, T_TP, DT_UP, DT_DN, DT_ST, DT_CH, DT_BM, DT_TP, T_P_GL, T_M_GL };
   void SetModel(FNC f, std::function<double(double, double, double)> model);

   // down, up, strange, charm, bottom, top
   std::array<std::function<double(double, double, double)>, 6> T
       = {zero_function, zero_function, zero_function, zero_function, zero_function, zero_function};

   std::array<std::function<double(double, double, double)>, 6> DT
       = {zero_function, zero_function, zero_function, zero_function, zero_function, zero_function};

   std::function<double(double, double, double)> T_p_gl = zero_function;
   std::function<double(double, double, double)> T_m_gl = zero_function;
}; // namespace Honeycomb

struct Solution {

   Solution(const Discretization *discretization, const InputModel &models, size_t nf);
   Solution() : _discretization(nullptr), _curr_basis(EVO)
   {
   }
   void PushFlavor();

   // Return the pairs of S^+_h, S^-_h, where h is popped flavor
   std::pair<Eigen::VectorXd, Eigen::VectorXd> PopFlavor();

   // Gl, S, NS1, ... => Gl, d, u, ...
   void RotateToPhysicalBasis();

   // Gl, d, u, ... => Gl, S, NS1, ...
   void RotateToEvolutionBasis();

   // Runge-Kutta methods
   void _copy(const Solution &other);
   void _plus_eq(double x, const Solution &other);
   void _ker_mul(double pref, const Kernels &ker);

public:
   const Discretization *_discretization;

   enum CURR_BASIS { PHYS, EVO };
   CURR_BASIS _curr_basis;

   size_t nf;

   // Gluon, Singlet, NS_1, NS_2, ...
   std::vector<Eigen::VectorXd> _distr_p;
   std::vector<Eigen::VectorXd> _distr_m;
};

// Thresholds in \mu^2
// returns vetor of intermediate scales between Q0 and Qf as
// {log(\mu_1^2), ..., \log(Qf^2)}. Q0 is not in the vector!
std::pair<std::vector<double>, Solution> get_initial_solution(double Q02, double Qf2,
                                                              const std::array<double, 6> &thresholds,
                                                              const Discretization *discretization,
                                                              const InputModel &models);

} // namespace Honeycomb

#endif // SOLUTION_HPP