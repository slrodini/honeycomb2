#include <honeycomb2/thread_pool.hpp>
#include <honeycomb2/utilities.hpp>
#include <honeycomb2/solution.hpp>
#include <honeycomb2/runge_kutta.hpp>
#include <honeycomb2/Eigen/Core>

namespace
{
unsigned int _private_num_av_threads = std::thread::hardware_concurrency();
unsigned int _private_num_threads    = _private_num_av_threads <= 2 ? 1 : _private_num_av_threads - 2;
Honeycomb::ThreadPool global_pool(_private_num_threads);
} // namespace

namespace Honeycomb
{

void InputModel::SetModel(FNC f, std::function<double(double, double, double)> model)
{
   switch (f) {
   case T_DN:
      T[0] = model;
      break;
   case T_UP:
      T[1] = model;
      break;
   case T_ST:
      T[2] = model;
      break;
   case T_CH:
      T[3] = model;
      break;
   case T_BM:
      T[4] = model;
      break;
   case T_TP:
      T[5] = model;
      break;
   case DT_DN:
      DT[0] = model;
      break;
   case DT_UP:
      DT[1] = model;
      break;
   case DT_ST:
      DT[2] = model;
      break;
   case DT_CH:
      DT[3] = model;
      break;
   case DT_BM:
      DT[4] = model;
      break;
   case DT_TP:
      DT[5] = model;
      break;
   case T_P_GL:
      T_p_gl = model;
      break;
   case T_M_GL:
      T_m_gl = model;
      break;
   default:
      break;
   }
}

Solution::Solution(const Discretization *discretization, const InputModel &models, size_t nf)
    : _discretization(discretization), nf(nf)
{

   std::vector<Eigen::VectorXd> v_Splus, v_Sminus;

   auto Fplus = [&models](double x1, double x2, double x3) -> double {
      return models.T_p_gl(x1, x2, x3) - models.T_p_gl(x1, x3, x2) + models.T_p_gl(x2, x1, x3);
   };

   auto Fminus = [&models](double x1, double x2, double x3) -> double {
      return models.T_m_gl(x1, x2, x3) + models.T_m_gl(x1, x3, x2) - models.T_m_gl(x2, x1, x3);
   };

   _distr_p.emplace_back(_discretization->discretize(Fplus));
   _distr_m.emplace_back(_discretization->discretize(Fminus));

   for (size_t i = 0; i < nf; i++) {
      auto Splus = [&models, i](double x1, double x2, double x3) -> double {
         return models.T[i](x1, x2, x3) - models.DT[i](x1, x2, x3) + models.T[i](x3, x2, x1)
              + models.DT[i](x3, x2, x1);
      };
      auto Sminus = [&models, i](double x1, double x2, double x3) -> double {
         return models.T[i](x1, x2, x3) - models.DT[i](x1, x2, x3) - models.T[i](x3, x2, x1)
              - models.DT[i](x3, x2, x1);
      };
      _distr_p.push_back(_discretization->discretize(Splus));
      _distr_m.push_back(_discretization->discretize(Sminus));
   }

   _curr_basis = PHYS;

   RotateToEvolutionBasis();
};

void Solution::PushFlavor()
{
   logger(Logger::INFO, "Pushing new falvor");
   // Trivial case
   if (nf == 1) {
      _distr_p.push_back(_distr_p[1]);
      _distr_m.push_back(_distr_m[1]);
      nf++;
      return;
   }
   // At the threshold, the new NS is equal to the current heaviest quark
   // for instance, at the charm threshold
   // NS3 = s - c, but c=0, so NS3 = s
   // Idem for all other cases.
   // To get the current heaviest flavor the general formula is
   // h = (Singlet - NS_{1} - 2 NS_{2} - ... - (nf-1)NS_{nf-1}) / nf

   // Get the singlet
   Eigen::VectorXd h_p = _distr_p[1];
   Eigen::VectorXd h_m = _distr_m[1];

   const double norm = 1.0 / static_cast<double>(nf);

   for (size_t i = 0; i < nf - 1; i++) {
      h_p -= (i + 1) * _distr_p[2 + i];
      h_m -= (i + 1) * _distr_m[2 + i];
   }
   _distr_p.push_back(h_p * norm);
   _distr_m.push_back(h_m * norm);

   nf++;
}

std::pair<Eigen::VectorXd, Eigen::VectorXd> Solution::PopFlavor()
{
   if (nf == 1) {
      logger(Logger::ERROR, "Solution::PopFlavor cannot remove flavor if only one is present.");
   }

   // I extract the flavor that will be killed and return it

   Eigen::VectorXd h_p = _distr_p[1];
   Eigen::VectorXd h_m = _distr_m[1];

   const double norm = 1.0 / static_cast<double>(nf);

   for (size_t i = 0; i < nf - 1; i++) {
      h_p -= (i + 1) * _distr_p[2 + i];
      h_m -= (i + 1) * _distr_m[2 + i];
   }

   _distr_p.pop_back();
   _distr_m.pop_back();

   // free unused memory
   _distr_p.shrink_to_fit();
   _distr_m.shrink_to_fit();

   nf--;

   return {h_p * norm, h_m * norm};
}

// Runge-Kutta methods
void Solution::_copy(const Solution &other)
{
   _distr_p    = other._distr_p;
   _distr_m    = other._distr_m;
   nf          = other.nf;
   _curr_basis = other._curr_basis;
}
void Solution::_plus_eq(double x, const Solution &other)
{
   // TODO: implement S += x*other
   size_t n = _distr_p.size();
   for (size_t i = 0; i < n; i++) {
      _distr_p[i] += x * other._distr_p[i];
      _distr_m[i] += x * other._distr_m[i];
   }
}
void Solution::_ker_mul(double pref, const Kernels &ker)
{
   // TODO: implement S = pref * (ker * S);
   global_pool.AddTask([&, this]() {
      for (size_t i = 2; i < _distr_p.size(); i++) {
         _distr_p[i] = pref * (ker.H_NS * _distr_p[i] - 3 * ker.CF * _distr_p[i]);
         _distr_m[i] = pref * (ker.H_NS * _distr_m[i] - 3 * ker.CF * _distr_m[i]);
      }
   });

   const double beta0 = (11.0 * ker.Nc - 2.0 * nf) / 3.0;

   global_pool.AddTask([&, this]() {
      Eigen::VectorXd tmp_g = _distr_p[0];
      Eigen::VectorXd tmp_s = _distr_p[1];

      _distr_p[0] = pref * (ker.H_gg_p * tmp_g - beta0 * tmp_g + ker.H_gq_p * tmp_s);
      _distr_p[1]
          = pref * (ker.H_NS * tmp_s - 3 * ker.CF * tmp_s + nf * ker.H_d13 * tmp_s + nf * ker.H_qg_p * tmp_g);

      tmp_g = _distr_m[0];
      tmp_s = _distr_m[1];

      _distr_m[0] = pref * (ker.H_gg_m * tmp_g - beta0 * tmp_g + ker.H_gq_m * tmp_s);
      _distr_m[1] = pref * (ker.H_NS * tmp_s - 3 * ker.CF * tmp_s + nf * ker.H_qg_m * tmp_g);
   });

   global_pool.WaitOnJobs();
}

std::pair<std::vector<double>, Solution> get_initial_solution(double Q02, double Qf2,
                                                              const std::array<double, 6> &thresholds,
                                                              const Discretization *discretization,
                                                              const InputModel &models)
{

   for (size_t i = 0; i < 5; i++) {
      if (thresholds[i + 1] < thresholds[i]) {
         logger(Logger::ERROR,
                "get_initial_solution: Thresholds are not increasing, backward evolution not supported yet.");
      }
   }
   const double tf = log(Qf2);

   std::vector<double> intermediate_scales;
   size_t nf = 0;
   for (size_t i = 0; i < thresholds.size(); i++) {
      if (thresholds[i] <= Q02) {
         nf++;
      } else if (thresholds[i] >= Qf2) {
         break;
      } else {
         intermediate_scales.push_back(log(thresholds[i]));
      }
   }
   intermediate_scales.push_back(tf);
   logger(Logger::INFO, std::format("nf at initial scale:      {:d}", nf));
   logger(Logger::INFO, std::format("# of intermediate scales: {:d}", intermediate_scales.size()));

   return {intermediate_scales, Solution(discretization, models, nf)};
}

void Solution::RotateToPhysicalBasis()
{
   if (_curr_basis == PHYS) {
      logger(Logger::WARNING, "Solution::RotateToPhysicalBasis already in physical basis. I do nothing.");
      return;
   }
   std::vector<Eigen::VectorXd> tmp_p;
   std::vector<Eigen::VectorXd> tmp_m;

   // Copy gluon
   tmp_p.push_back(_distr_p[0]);
   tmp_m.push_back(_distr_m[0]);

   if (nf != _distr_p.size() - 1) {
      logger(Logger::ERROR, "Bug in push flavor, vector of solution has incorrect dimensions.");
   }
   int nfi = static_cast<int>(nf);

   // Construct from lightest to heaviest
   for (int k = 0; k < nfi; k++) {
      Eigen::VectorXd h_p = _distr_p[1]; // start with h = S
      Eigen::VectorXd h_m = _distr_m[1]; // start with h = S
      for (int i = nfi - 2; i >= k; i--) {
         h_p += (nfi - 1 - i) * _distr_p[2 + i];
         h_m += (nfi - 1 - i) * _distr_m[2 + i];
      }
      for (int i = 0; i < k; i++) {
         h_p -= (i + 1) * _distr_p[2 + i];
         h_m -= (i + 1) * _distr_m[2 + i];
      }
      tmp_p.push_back(h_p / nfi);
      tmp_m.push_back(h_m / nfi);
   }
   _distr_p    = tmp_p;
   _distr_m    = tmp_m;
   _curr_basis = PHYS;
}

void Solution::RotateToEvolutionBasis()
{
   if (_curr_basis == EVO) {
      logger(Logger::WARNING, "Solution::RotateToEvolutionBasis already in evolution basis. I do nothing.");
      return;
   }
   std::vector<Eigen::VectorXd> tmp_p;
   std::vector<Eigen::VectorXd> tmp_m;

   // Copy gluon
   tmp_p.push_back(_distr_p[0]);
   tmp_m.push_back(_distr_m[0]);

   int nfi = static_cast<int>(nf);

   // Construct from lightest to heaviest
   Eigen::VectorXd h_p = _distr_p[1]; // start with h = d
   Eigen::VectorXd h_m = _distr_m[1]; // start with h = d
   for (int k = 1; k < nfi; k++) {
      h_p += _distr_p[1 + k];
      h_m += _distr_m[1 + k];
   }
   tmp_p.push_back(h_p);
   tmp_m.push_back(h_m);

   for (int k = 0; k < nfi - 1; k++) {
      tmp_p.push_back(_distr_p[1 + k] - _distr_p[2 + k]);
      tmp_m.push_back(_distr_m[1 + k] - _distr_m[2 + k]);
   }
   _distr_p    = tmp_p;
   _distr_m    = tmp_m;
   _curr_basis = EVO;
}

EvolutionOperatorFixedNf::EvolutionOperatorFixedNf(const Grid2D *grid, size_t _nf) : _grid(grid), nf(_nf)
{

   unsigned int num_av_threads = std::thread::hardware_concurrency();
   unsigned int num_threads    = num_av_threads <= 2 ? 1 : num_av_threads - 2;

   th_pool = new ThreadPool(num_threads);

   NS_P = Eigen::MatrixXd::Identity(grid->c_size_li, grid->c_size_li);
   NS_M = Eigen::MatrixXd::Identity(grid->c_size_li, grid->c_size_li);

   S_P = Eigen::MatrixXd::Identity(2 * grid->c_size_li, 2 * grid->c_size_li);
   S_M = Eigen::MatrixXd::Identity(2 * grid->c_size_li, 2 * grid->c_size_li);
}

void EvolutionOperatorFixedNf::_copy(const EvolutionOperatorFixedNf &other)
{
   NS_P = other.NS_P;
   NS_M = other.NS_M;
   S_P  = other.S_P;
   S_M  = other.S_M;

   th_pool = other.th_pool;
}

void EvolutionOperatorFixedNf::_plus_eq(double x, const EvolutionOperatorFixedNf &other)
{
   NS_P += x * other.NS_P;
   NS_M += x * other.NS_M;
   S_P  += x * other.S_P;
   S_M  += x * other.S_M;
}

void EvolutionOperatorFixedNf::_ker_mul(double pref, const MergedKernelsFixedNf &ker)
{
   th_pool->AddTask([&]() {
      NS_P = pref * (ker.H_NS * NS_P);
   });
   th_pool->AddTask([&]() {
      NS_M = pref * (ker.H_NS * NS_M);
   });
   th_pool->AddTask([&]() {
      S_P = pref * (ker.H_S_P * S_P);
   });
   th_pool->AddTask([&]() {
      S_M = pref * (ker.H_S_M * S_M);
   });
   th_pool->WaitOnJobs();
}

void ApplyEvolutionOperator(Solution &sol, const EvolutionOperatorFixedNf &O)
{
   if (sol.nf != O.nf) {
      logger(Logger::ERROR, std::format("ApplyEvolutionOperator: incompatible nf between Solution ({:d}) and "
                                        "Evolution Operator ({:d}).",
                                        sol.nf, O.nf));
   }

   for (size_t i = 2; i < sol._distr_p.size(); i++) {
      sol._distr_p[i] = O.NS_P * sol._distr_p[i];
      sol._distr_m[i] = O.NS_M * sol._distr_m[i];
   }

   const long int r      = sol._distr_p[0].size();
   Eigen::VectorXd tmp_p = Eigen::VectorXd::Zero(r * 2);
   Eigen::VectorXd tmp_m = Eigen::VectorXd::Zero(r * 2);

   for (long int i = 0; i < r; i++) {
      tmp_p(i)     = sol._distr_p[0](i);
      tmp_p(i + r) = sol._distr_p[1](i);

      tmp_m(i)     = sol._distr_m[0](i);
      tmp_m(i + r) = sol._distr_m[1](i);
   }

   tmp_p = O.S_P * tmp_p;
   tmp_m = O.S_M * tmp_m;

   for (long int i = 0; i < r; i++) {
      sol._distr_p[0](i) = tmp_p(i);
      sol._distr_p[1](i) = tmp_p(i + r);
      sol._distr_m[0](i) = tmp_m(i);
      sol._distr_m[1](i) = tmp_m(i + r);
   }
}

Solution evolve_solution(const Kernels &kers, double Q02, double Qf2, const std::array<double, 6> &thresholds,
                         const Discretization *discretization, const InputModel &models,
                         std::function<double(double)> as)
{
   auto [inter_scales, sol0] = get_initial_solution(Q02, Qf2, thresholds, discretization, models);

   const long int r = sol0._distr_p[0].size();

   inter_scales.insert(inter_scales.begin(), log(Q02));

   size_t nf_fin = 0;
   size_t nf_in  = 0;
   for (size_t i = 0; i < thresholds.size(); i++) {
      if (thresholds[i] <= Qf2) {
         nf_fin++;
      }

      if (thresholds[i] <= Q02) {
         nf_in++;
      }
   }

   logger(Logger::INFO, std::format("nf in {:d}, nf end {:d}", nf_in, nf_fin));

   std::vector<MergedKernelsFixedNf> nf_kers;
   for (size_t nf = nf_in; nf <= nf_fin; nf++) {
      nf_kers.emplace_back(MergedKernelsFixedNf(kers, nf));
   }
   std::vector<size_t> n_steps = {20};

   for (size_t i = 0; i < inter_scales.size() - 1; i++) {

      // Non-singlets
      for (size_t j = 2; j < sol0._distr_p.size(); j++) {
         global_pool.AddTask([&, j]() {
            runge_kutta::GenericRungeKutta<Eigen::MatrixXd, Eigen::VectorXd, 13> evolver_p(
                nf_kers[i].H_NS, sol0._distr_p[j], Honeycomb::runge_kutta::DOPRI8, as, -1.0, inter_scales[i],
                0.01);
            evolver_p({inter_scales[i + 1]}, n_steps);
            sol0._distr_p[j] = evolver_p.GetSolution();
         });

         global_pool.AddTask([&, j]() {
            runge_kutta::GenericRungeKutta<Eigen::MatrixXd, Eigen::VectorXd, 13> evolver_m(
                nf_kers[i].H_NS, sol0._distr_m[j], Honeycomb::runge_kutta::DOPRI8, as, -1.0, inter_scales[i],
                0.01);
            evolver_m({inter_scales[i + 1]}, n_steps);
            sol0._distr_m[j] = evolver_m.GetSolution();
         });
      }

      global_pool.AddTask([&]() {
         Eigen::VectorXd tmp_p = Eigen::VectorXd::Zero(r * 2);
         Eigen::VectorXd tmp_m = Eigen::VectorXd::Zero(r * 2);

         for (long int k = 0; k < r; k++) {
            tmp_p(k)     = sol0._distr_p[0](k);
            tmp_p(k + r) = sol0._distr_p[1](k);

            tmp_m(k)     = sol0._distr_m[0](k);
            tmp_m(k + r) = sol0._distr_m[1](k);
         }

         runge_kutta::GenericRungeKutta<Eigen::MatrixXd, Eigen::VectorXd, 13> evolver_p(
             nf_kers[i].H_S_P, tmp_p, Honeycomb::runge_kutta::DOPRI8, as, -1.0, inter_scales[i], 0.01);
         evolver_p({inter_scales[i + 1]}, n_steps);
         tmp_p = evolver_p.GetSolution();

         runge_kutta::GenericRungeKutta<Eigen::MatrixXd, Eigen::VectorXd, 13> evolver_m(
             nf_kers[i].H_S_M, tmp_m, Honeycomb::runge_kutta::DOPRI8, as, -1.0, inter_scales[i], 0.01);
         evolver_m({inter_scales[i + 1]}, n_steps);
         tmp_m = evolver_m.GetSolution();

         for (long int k = 0; k < r; k++) {
            sol0._distr_p[0](k) = tmp_p(k);
            sol0._distr_p[1](k) = tmp_p(k + r);
            sol0._distr_m[0](k) = tmp_m(k);
            sol0._distr_m[1](k) = tmp_m(k + r);
         }
      });
      global_pool.WaitOnJobs();

      if (i != inter_scales.size() - 2) sol0.PushFlavor();
   }

   return sol0;
}

} // namespace Honeycomb