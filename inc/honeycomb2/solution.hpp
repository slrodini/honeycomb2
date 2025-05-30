#ifndef HC2_SOLUTION_HPP
#define HC2_SOLUTION_HPP

#include <honeycomb2/default.hpp>
#include <honeycomb2/kernels.hpp>
#include <honeycomb2/discretization.hpp>
#include <honeycomb2/random_engine.hpp>
#include <honeycomb2/thread_pool.hpp>
#include <honeycomb2/Eigen/Core>

namespace Honeycomb
{

inline double zero_function(double, double, double)
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

   // Utilities
   bool is_equalt_to(const Solution &other, double acc = 1.0e-12) const;

public:
   const Discretization *_discretization;

   enum CURR_BASIS { PHYS, EVO };
   CURR_BASIS _curr_basis;

   size_t nf;

   // Gluon, Singlet, NS_1, NS_2, ...
   std::vector<Eigen::VectorXd> _distr_p;
   std::vector<Eigen::VectorXd> _distr_m;
};

struct OutputModel {
   enum FNC { T_UP, T_DN, T_ST, T_CH, T_BM, T_TP, DT_UP, DT_DN, DT_ST, DT_CH, DT_BM, DT_TP, T_P_GL, T_M_GL };

   OutputModel(const Solution &sol);

   double GetDistribution(FNC f, const RnC::Triplet &x123)
   {
      return GetDistribution(f, RnC::from_x123_to_rhophi(x123));
   }
   double GetDistribution(FNC f, const double x1, const double x2, const double x3)
   {
      return GetDistribution(f, RnC::from_x123_to_rhophi({x1, x2, x3}));
   }
   double GetDistribution(FNC f, const RnC::Pair &rhophi);

   // Note: First element is gluon T^+ and T^-, not DeltaT!
   std::vector<Eigen::VectorXd> T;
   std::vector<Eigen::VectorXd> DT;
   const Discretization *_discretization;
};

// Struct to be fed to Runge-Kutta
struct EvolutionOperatorFixedNf {

   EvolutionOperatorFixedNf(const Grid2D *grid, size_t nf, ThreadPool *th);
   EvolutionOperatorFixedNf() : _grid(nullptr), th_pool(nullptr)
   {
   }

   void PushFlavor();

   // Runge-Kutta methods
   void _copy(const EvolutionOperatorFixedNf &other);
   void _plus_eq(double x, const EvolutionOperatorFixedNf &other);
   void _ker_mul(double pref, const MergedKernelsFixedNf &ker);

public:
   const Grid2D *_grid;
   ThreadPool *th_pool;

   size_t nf;

   Eigen::MatrixXd NS_P, S_P;
   Eigen::MatrixXd NS_M, S_M;
};

// Final struct, simplified structure
struct EvOpNF {
   EvOpNF(const EvolutionOperatorFixedNf &);
   EvOpNF() = default;

   void SetScales(double t0, double tF);

   size_t nf;
   Eigen::MatrixXd NS_P, S_P;
   Eigen::MatrixXd NS_M, S_M;
   std::pair<double, double> t0tF;

   template <class Archive>
   void save(Archive &archive) const
   {
      archive(nf);

      archive(t0tF);

      // _P and _M same dimensions
      archive(NS_P.rows());
      archive(NS_P.cols());

      archive(S_P.rows());
      archive(S_P.cols());

      // Explicitly archive matrix elements
      for (int i = 0; i < NS_P.rows(); ++i) {
         for (int j = 0; j < NS_P.cols(); ++j) {
            archive(NS_P(i, j));
            archive(NS_M(i, j));
         }
      }
      for (int i = 0; i < S_P.rows(); ++i) {
         for (int j = 0; j < S_P.cols(); ++j) {
            archive(S_P(i, j));
            archive(S_M(i, j));
         }
      }
   }

   template <class Archive>
   void load(Archive &archive)
   {
      archive(nf);

      archive(t0tF);

      // _P and _M same dimensions
      long int ns_r, ns_c, s_r, s_c;
      archive(ns_r);
      archive(ns_c);

      archive(s_r);
      archive(s_c);
      NS_P = Eigen::MatrixXd::Zero(ns_r, ns_c);
      NS_M = Eigen::MatrixXd::Zero(ns_r, ns_c);

      S_P = Eigen::MatrixXd::Zero(s_r, s_c);
      S_M = Eigen::MatrixXd::Zero(s_r, s_c);
      for (int i = 0; i < NS_P.rows(); ++i) {
         for (int j = 0; j < NS_P.cols(); ++j) {
            archive(NS_P(i, j));
            archive(NS_M(i, j));
         }
      }
      for (int i = 0; i < S_P.rows(); ++i) {
         for (int j = 0; j < S_P.cols(); ++j) {
            archive(S_P(i, j));
            archive(S_M(i, j));
         }
      }
   }
};

struct EvOp {
   EvOp(Grid2D *grid) : _grid(grid) {};
   EvOp() : _grid(nullptr) {};

   void emplace_back(const EvOpNF &arg);
   void push_back(const EvOpNF &arg);
   void SetScales(double t0, double tF);
   void SetThresholds(const std::vector<double> &thresholds);

   EvOpNF &back()
   {
      return _operators.back();
   }

   Grid2D *_grid;
   std::vector<EvOpNF> _operators;
   std::vector<double> _thresholds; // This stores ALL the 6 thresholds
   std::pair<double, double> t0tF;
   template <class Archive>
   void save(Archive &archive) const
   {
      archive(*_grid);
      archive(_operators);
      archive(_thresholds);
      archive(t0tF);
   }

   template <class Archive>
   void load(Archive &archive)
   {
      Grid2D tmp_grid;
      archive(tmp_grid);
      check_grid_compatibility(tmp_grid, *_grid);

      archive(_operators);
      archive(_thresholds);
      archive(t0tF);
   }
};

inline void save_evolution_operator(const EvOp &O, const std::string &file_name)
{
   SaveChecksumArchive<EvOp, cereal::PortableBinaryOutputArchive>(O, file_name);
}

inline std::pair<bool, EvOp> load_evolution_operator(const std::string &file_name, Grid2D *g)
{
   EvOp result(g);
   if (!LoadAndVerify<EvOp, cereal::PortableBinaryInputArchive>(file_name, result)) {
      logger(Logger::WARNING, "I was not able to correctly load the cereal archive " + file_name
                                  + " containing the evolution operator.");
      return {false, result};
   }
   return {true, result};
}

EvOp compute_evolution_operator(Grid2D *grid, const Kernels &kers, double Q02, double Qf2,
                                const std::array<double, 6> &thresholds, std::function<double(double)> as);

// Thresholds in \mu^2
// returns vetor of intermediate scales between Q0 and Qf as
// {log(\mu_1^2), ..., \log(Qf^2)}. Q0 is not in the vector!
std::pair<std::vector<double>, Solution> get_initial_solution(double Q02, double Qf2,
                                                              const std::array<double, 6> &thresholds,
                                                              const Discretization *discretization,
                                                              const InputModel &models);

void ApplyEvolutionOperator(Solution &sol, const EvOpNF &O);
void ApplyEvolutionOperator(Solution &sol, const EvOp &O);

Solution evolve_solution(const Kernels &kers, double Q02, double Qf2, const std::array<double, 6> &thresholds,
                         const Discretization *discretization, const InputModel &models,
                         std::function<double(double)> as);

} // namespace Honeycomb

#endif // HC2_SOLUTION_HPP