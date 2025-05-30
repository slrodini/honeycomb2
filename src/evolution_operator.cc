#include <honeycomb2/discretization.hpp>
#include <honeycomb2/thread_pool.hpp>
#include <honeycomb2/utilities.hpp>
#include <honeycomb2/solution.hpp>
#include <honeycomb2/runge_kutta.hpp>
#include <honeycomb2/Eigen/Core>

namespace Honeycomb
{
EvOpNF::EvOpNF(const EvolutionOperatorFixedNf &O)
    : nf(O.nf), NS_P(O.NS_P), S_P(O.S_P), NS_M(O.NS_M), S_M(O.S_M)
{
}

void EvOpNF::SetScales(double t0, double tF)
{
   t0tF.first  = t0;
   t0tF.second = tF;
}

void EvOp::SetScales(double t0, double tF)
{
   t0tF.first  = t0;
   t0tF.second = tF;
}

void EvOp::SetThresholds(const std::vector<double> &thresholds)
{
   _thresholds = thresholds;
}

void EvOp::emplace_back(const EvOpNF &arg)
{
   _operators.emplace_back(arg);
}

void EvOp::push_back(const EvOpNF &arg)
{
   _operators.push_back(arg);
}

// -----------------------------------------------------------------------------

EvOp compute_evolution_operator(Grid2D *grid, const Kernels &kers, double Q02, double Qf2,
                                const std::array<double, 6> &thresholds, std::function<double(double)> as)
{

   EvOp result(grid);
   result.SetScales(log(Q02), log(Qf2));
   result.SetThresholds(std::vector<double>(thresholds.begin(), thresholds.end()));

   std::vector<double> inter_scales;

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

   for (size_t nf = nf_in; nf < nf_fin; nf++) {
      inter_scales.push_back(log(thresholds[nf]));
   }
   inter_scales.insert(inter_scales.begin(), log(Q02));
   inter_scales.push_back(log(Qf2));
   // inter_scales = [t_{Q0}, t_{mq1}, t_{mq2}, ..., t_{Qf}]
   // t_x = log(x**2)

   logger(Logger::INFO, std::format("nf in {:d}, nf end {:d}", nf_in, nf_fin));

   std::vector<MergedKernelsFixedNf> nf_kers;
   for (size_t nf = nf_in; nf <= nf_fin; nf++) {
      nf_kers.emplace_back(MergedKernelsFixedNf(kers, nf));
   }
   std::vector<size_t> n_steps = {20};
   unsigned int num_av_threads = std::thread::hardware_concurrency();
   unsigned int num_threads    = num_av_threads <= 2 ? 1 : num_av_threads - 2;
   ThreadPool local_pool(num_threads);

   size_t curr_nf = nf_in;
   for (size_t i = 0; i < inter_scales.size() - 1; i++, curr_nf++) {

      EvolutionOperatorFixedNf O1(&kers.grid, curr_nf, &local_pool);

      Honeycomb::runge_kutta::GenericRungeKutta<Honeycomb::MergedKernelsFixedNf,
                                                Honeycomb::EvolutionOperatorFixedNf, 13>
          evolver_O(nf_kers[curr_nf - nf_in], O1, Honeycomb::runge_kutta::DOPRI8, as, -1.0, inter_scales[i],
                    0.01);

      evolver_O({inter_scales[i + 1]}, n_steps);
      result.push_back(EvOpNF(evolver_O.GetSolution()));
      result.back().SetScales(inter_scales[i], inter_scales[i + 1]);
   }

   return result;
}

void ApplyEvolutionOperator(Solution &sol, const EvOp &O)
{
   for (size_t i = 0; i < O._operators.size(); i++) {
      ApplyEvolutionOperator(sol, O._operators[i]);

      if (i != O._operators.size() - 1) sol.PushFlavor();
   }
}

} // namespace Honeycomb