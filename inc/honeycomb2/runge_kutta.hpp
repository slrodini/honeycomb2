#ifndef HC2_RUNGE_KUTTA_HPP
#define HC2_RUNGE_KUTTA_HPP

#include <functional>
#include <concepts>
#include <type_traits>
#include <cmath>

#include <iostream>
#include <vector>
#include <array>

#include <honeycomb2/utilities.hpp>

namespace Honeycomb
{
namespace rk
{
/**
 * @brief Butcher Tableau for a generic Runge-Kutta method
 *
 * @tparam Stages
 */
template <size_t Stages>
struct ButcherTableau {

   /// Delete default constructor to avoid non-sense
   ButcherTableau() = delete;

   ButcherTableau(double order, const std::array<double, Stages> &ci,
                  const std::array<double, Stages> &bi,
                  const std::array<std::vector<double>, Stages> &ai)
       : _order(order), _aij(ai), _bi(bi), _ci(ci)
   {
      double s = 0.;
      /// Sanity checks
      for (size_t i = 0; i < Stages; i++) {
         s += bi[i];
         if (ai[i].size() > i)
            throw std::invalid_argument(
                "Implicit methods (for which aij != 0  if j>=i) are not supported.");
      }
      if (std::fabs(s - 1.) > 1.0e-14) {
         std::printf("%.16e\n", s);
         throw std::invalid_argument("ButcherTableau: Sum of bi must be 1;");
      }
   };

   static constexpr size_t stages = Stages;
   const double _order;
   const std::array<std::vector<double>, Stages> _aij;
   const std::array<double, Stages> _bi;
   const std::array<double, Stages> _ci;
};

template <size_t Stages>
struct ButcherTableauWErrorEstimate : public ButcherTableau<Stages> {

   /// Delete default constructor to avoid non-sense
   ButcherTableauWErrorEstimate() = delete;

   ButcherTableauWErrorEstimate(double hi, double lo, const std::array<double, Stages> &ci,
                                const std::array<double, Stages> &bi,
                                const std::array<double, Stages> &bi_lo,
                                const std::array<std::vector<double>, Stages> &ai)
       : ButcherTableau<Stages>(hi, ci, bi, ai), _order_lo(lo), _bi_lo(bi_lo) {
            // Note: depending on the method, not all bi_lo will sum to 1.
         };

   const double _order_lo;
   const std::array<double, Stages> _bi_lo;
};

/// Minimal required interface
template <typename S>
concept RungeKuttaCompatible = requires {
   // Exact match for: void clone(const S& other)
   { static_cast<void (S::*)(const S &)>(&S::clone) };

   // Exact match for: void add_with_weight(double, const S&)
   // Must implement S += w * S1, i.e. S = S + w * S1
   { static_cast<void (S::*)(double, const S &)>(&S::add_with_weight) };

   // Multiply by scalar
   { static_cast<void (S::*)(double)>(&S::scalar_mult) };
};

struct TimeInfo {
   TimeInfo(const std::vector<double> &ts = {0., 1.}, double dt_min = 0.1, size_t n_step_min = 10)
   {
      if (ts.size() <= 1) {
         throw std::invalid_argument(
             "TimeInfo: at least 2 time points are required (initial and final)");
      }
      _ts = ts;
      std::sort(_ts.begin(), _ts.end());

      double nd = static_cast<double>(n_step_min);
      for (size_t i = 0; i < _ts.size() - 1; i++) {
         size_t n_step = n_step_min;

         double dt_tmp = (_ts[i + 1] - _ts[i]) / nd;
         while (dt_tmp > dt_min) {
            n_step *= 2;
            dt_tmp  = (_ts[i + 1] - _ts[i]) / static_cast<double>(n_step);
         }
         _dt.push_back(dt_tmp);
         _n_step.push_back(n_step);
      }
   }

   std::vector<double> _ts, _dt;
   std::vector<size_t> _n_step;
};

/// Useful aliases
template <typename Solution>
using rk_rhs_t = std::function<void(double, Solution &)>;

template <typename Solution>
using rk_callback_t = std::function<void(double, const TimeInfo &, Solution &)>;

template <typename Solution, size_t Stages>
requires RungeKuttaCompatible<Solution>
struct RungeKutta {
public:
   RungeKutta(const ButcherTableau<Stages> &tableau, const Solution &initial_conditions,
              TimeInfo time_info, rk_rhs_t<Solution> &rhs)
       : _cb_each_step(false), _tableau(tableau), _solution(initial_conditions),
         _temp_step(initial_conditions), _time_info(time_info), _rhs(rhs)
   {
      for (size_t i = 0; i < Stages; i++) {
         _ki.push_back(initial_conditions);
      }
   }

   void AddCallback(rk_callback_t<Solution> &callback)
   {
      _callback = callback;
   }

   void CallbackEachStep()
   {
      _cb_each_step = true;
   }

   void CallbackOnlyOnTimeStamp()
   {
      _cb_each_step = false;
   }

   void step()
   {

      _temp_step.clone(_solution);

      for (size_t i = 0; i < Stages; i++) {

         _ki[i].clone(_solution);

         for (size_t j = 0; j < i; j++) {
            if (_tableau._aij[i][j] == 0.) continue;
            _ki[i].add_with_weight(_tableau._aij[i][j], _ki[j]);
         }

         _rhs(_t + _dt * _tableau._ci[i], _ki[i]);
         _ki[i].scalar_mult(_dt);

         _temp_step.add_with_weight(_tableau._bi[i], _ki[i]);
      }

      _solution.clone(_temp_step);
   }

   /// @name For visualization applications
   /// @{
   void advance_t()
   {
      _t += _dt;
   }
   void set_dt(double dt)
   {
      _dt = dt;
   }
   /// @}

   void operator()()
   {
      for (size_t i = 0; i < _time_info._ts.size() - 1; i++) {
         _t = _time_info._ts[i];
         set_dt(_time_info._dt[i]);
         for (size_t j = 0; j < _time_info._n_step[i]; j++) {
            step();
            advance_t();
            if (_callback && _cb_each_step) {
               (*_callback)(_t, _time_info, _solution);
            }
         }

         if (_callback && !_cb_each_step) {
            (*_callback)(_t, _time_info, _solution);
         }

         if (std::fabs(_t - _time_info._ts[i + 1]) > 1.0e-12) {
            std::cerr << std::fabs(_t - _time_info._ts[i + 1]) << std::endl;
            throw std::runtime_error("RungeKutta: time steps do not match. This is a bug.");
         }
      }
   }

   const Solution &GetSolution()
   {
      return _solution;
   }

private:
   double _t, _dt;

   bool _cb_each_step;
   ButcherTableau<Stages> _tableau;
   Solution _solution;
   Solution _temp_step;
   std::vector<Solution> _ki;

   TimeInfo _time_info;

   /// The right hand side of ODE
   std::reference_wrapper<rk_rhs_t<Solution>> _rhs;
   /// The (optional) callback
   std::optional<std::reference_wrapper<rk_callback_t<Solution>>> _callback;
};

struct PreImplementedTableau {
   static inline const ButcherTableau<4> RKOriginal
       = ButcherTableau<4>(4, {0, 0.5, 0.5, 1},                  /* ci  */
                           {1. / 6., 1. / 3., 1. / 3., 1. / 6.}, /* bi  */
                           {{                                    /* aij */
                             {},
                             {0.5},
                             {0, 0.5},
                             {0, 0, 1}}});

   static inline const ButcherTableau<13> DOPRI8 = ButcherTableau<13>(
       8,
       {0., 1. / 18, 1. / 12, 1. / 8, 5. / 16, 3. / 8, 59. / 400, 93. / 200,
        5490023248. / 9719169821, 13. / 20, 1201146811. / 1299019798, 1., 1.},
       {14005451. / 335480064, 0., 0., 0., 0., -59238493. / 1068277825, 181606767. / 758867731,
        561292985. / 797845732, -1041891430. / 1371343529, 760417239. / 1151165299,
        118820643. / 751138087, -528747749. / 2220607170, 1. / 4},
       {{{},
         {1. / 18},
         {1. / 48, 1. / 16},
         {1. / 32, 0., 3. / 32},
         {5. / 16, 0., -75. / 64, 75. / 64},
         {3. / 80, 0., 0., 3. / 16, 3. / 20},
         {29443841. / 614563906, 0., 0., 77736538. / 692538347, -28693883. / 1125000000,
          23124283. / 1800000000},
         {16016141. / 946692911, 0., 0., 61564180. / 158732637, 22789713. / 633445777,
          545815736. / 2771057229, -180193667. / 1043307555},
         {39632708. / 573591083, 0., 0., -433636366. / 683701615, -421739975. / 2616292301,
          100302831. / 723423059, 790204164. / 839813087, 800635310. / 3783071287},
         {246121993. / 1340847787, 0., 0., -37695042795. / 15268766246, -309121744. / 1061227803,
          -12992083. / 490766935, 6005943493. / 2108947869, 393006217. / 1396673457,
          123872331. / 1001029789},
         {-1028468189. / 846180014, 0., 0., 8478235783. / 508512852, 1311729495. / 1432422823,
          -10304129995. / 1701304382, -48777925059. / 3047939560, 15336726248. / 1032824649,
          -45442868181. / 3398467696, 3065993473. / 597172653},
         {185892177. / 718116043, 0., 0., -3185094517. / 667107341, -477755414. / 1098053517,
          -703635378. / 230739211, 5731566787. / 1027545527, 5232866602. / 850066563,
          -4093664535. / 808688257, 3962137247. / 1805957418, 65686358. / 487910083},
         {403863854. / 491063109, 0., 0., -5068492393. / 434740067, -411421997. / 543043805,
          652783627. / 914296604, 11173962825. / 925320556, -13158990841. / 6184727034,
          3936647629. / 1978049680, -160528059. / 685178525, 248638103. / 1413531060, 0.}}});

   static inline const ButcherTableauWErrorEstimate<7> TSIT54 = ButcherTableauWErrorEstimate<7>(
       5, 4, {0., 0.161, 0.327, 0.9, 0.9800255409045097, 1.0, 1.0},
       {0.09646076681806523, 0.01, 0.4798896504144996, 1.379008574103742, -3.290069515436081,
        2.324710524099774, 0.},
       {0.001780011052226, 0.000816434459657, -0.007880878010262, 0.144711007173263,
        -0.582357165452555, 0.458082105929187, 1. / 66.},
       {{
           {},
           {0.161},
           {-0.008480655492356992, 0.3354806554923570},
           {2.8971530571054944, -6.359448489975075, 4.362295432869581},
           {5.32586482843926, -11.74888356406283, 7.495539342889836, -0.09249506636175525},
           {5.661747395595482, -12.92096931784711, 8.159367898576159, 0.07158497328140100,
            0.02826905039406838},
           {0.09646076681806523, 0.01, 0.4798896504144996, 1.379008574103742, -3.290069515436081,
            2.324710524099774},
       }});

   static inline const ButcherTableauWErrorEstimate<13> DOPRI87 = ButcherTableauWErrorEstimate<13>(
       8, 7,
       {0., 1. / 18, 1. / 12, 1. / 8, 5. / 16, 3. / 8, 59. / 400, 93. / 200,
        5490023248. / 9719169821, 13. / 20, 1201146811. / 1299019798, 1., 1.},
       {14005451. / 335480064, 0., 0., 0., 0., -59238493. / 1068277825, 181606767. / 758867731,
        561292985. / 797845732, -1041891430. / 1371343529, 760417239. / 1151165299,
        118820643. / 751138087, -528747749. / 2220607170, 1. / 4},
       {13451932. / 455176623., 0., 0., 0., 0., -808719846. / 976000145., 1757004468. / 5645159321.,
        656045339. / 265891186., -3867574721. / 1518517206., 465885868. / 322736535.,
        53011238. / 667516719., 2. / 45., 0.},
       {{{},
         {1. / 18},
         {1. / 48, 1. / 16},
         {1. / 32, 0., 3. / 32},
         {5. / 16, 0., -75. / 64, 75. / 64},
         {3. / 80, 0., 0., 3. / 16, 3. / 20},
         {29443841. / 614563906, 0., 0., 77736538. / 692538347, -28693883. / 1125000000,
          23124283. / 1800000000},
         {16016141. / 946692911, 0., 0., 61564180. / 158732637, 22789713. / 633445777,
          545815736. / 2771057229, -180193667. / 1043307555},
         {39632708. / 573591083, 0., 0., -433636366. / 683701615, -421739975. / 2616292301,
          100302831. / 723423059, 790204164. / 839813087, 800635310. / 3783071287},
         {246121993. / 1340847787, 0., 0., -37695042795. / 15268766246, -309121744. / 1061227803,
          -12992083. / 490766935, 6005943493. / 2108947869, 393006217. / 1396673457,
          123872331. / 1001029789},
         {-1028468189. / 846180014, 0., 0., 8478235783. / 508512852, 1311729495. / 1432422823,
          -10304129995. / 1701304382, -48777925059. / 3047939560, 15336726248. / 1032824649,
          -45442868181. / 3398467696, 3065993473. / 597172653},
         {185892177. / 718116043, 0., 0., -3185094517. / 667107341, -477755414. / 1098053517,
          -703635378. / 230739211, 5731566787. / 1027545527, 5232866602. / 850066563,
          -4093664535. / 808688257, 3962137247. / 1805957418, 65686358. / 487910083},
         {403863854. / 491063109, 0., 0., -5068492393. / 434740067, -411421997. / 543043805,
          652783627. / 914296604, 11173962825. / 925320556, -13158990841. / 6184727034,
          3936647629. / 1978049680, -160528059. / 685178525, 248638103. / 1413531060, 0.}}});

   /// Next three tableau from https://en.wikipedia.org/wiki/Runge-Kutta-Fehlberg_method
   static inline const ButcherTableauWErrorEstimate<6> FEHL_F2 = ButcherTableauWErrorEstimate<6>(
       5, 4, {0., 0.25, 3. / 8., 12. / 13., 1., 0.5},
       {16. / 135., 0., 6656. / 12825., 28561. / 56430., -9. / 50., 2. / 55.},
       {25. / 216., 0., 1408. / 2565., 2197 / 4104., -1. / 5., 0.},
       {{
           {},
           {1. / 4.},
           {3. / 32., 9. / 32.},
           {1932. / 2197., -7200. / 2197., 7296. / 2197.},
           {439. / 216., -8., 3680. / 513., -845. / 4104.},
           {-8. / 27., 2., -3544. / 2565., 1859. / 4104., -11. / 40.},
       }});

   static inline const ButcherTableauWErrorEstimate<6> FEHL_F1 = ButcherTableauWErrorEstimate<6>(
       5, 4, {0., 2. / 9., 1. / 3., 3. / 4., 1., 5. / 6.},
       {47. / 450., 0., 12. / 25., 32. / 225., 1. / 30., 6. / 25.},
       {1. / 9., 0., 9. / 20., 16. / 45., 1. / 12., 0.},
       {{
           {},
           {2. / 9.},
           {1. / 12., 1. / 4.},
           {69. / 128., -243. / 128., 135. / 64.},
           {-17. / 12., 27. / 4., -27. / 5., 16. / 15.},
           {65 / 432., -5. / 16., 13. / 16., 4. / 27., 5. / 144.},
       }});

   static inline const ButcherTableauWErrorEstimate<6> FEHL_F4 = ButcherTableauWErrorEstimate<6>(
       5, 4, {0., 0.5, 0.5, 1., 2. / 3., 1. / 5.},
       {1. / 24., 0., 0., 5. / 48., 27. / 56., 125. / 336.},
       {1. / 6., 0., 2. / 3., 1. / 6., 0., 0.},
       {{
           {},
           {0.5},
           {0.25, 0.25},
           {0., -1., 2.},
           {7. / 27., 10. / 27., 0., 1. / 27.},
           {28. / 625., -1. / 5., 546. / 625., 54. / 625., -378 / 625.},
       }});

   static inline const ButcherTableauWErrorEstimate<7> DOPRI54 = ButcherTableauWErrorEstimate<7>(
       5, 4, {0., 1. / 5, 3. / 10, 4. / 5, 8. / 9, 1., 1.},
       {35. / 384, 0., 500. / 1113, 125. / 192, -2187. / 6784, 11. / 84, 0.},
       {5179. / 57600., 0., 7571. / 16695., 393. / 640., -92097. / 339200., 187. / 2100., 1. / 40.},
       {{{},
         {1. / 5},
         {3. / 40, 9. / 40},
         {44. / 45, -56. / 15, 32. / 9},
         {19372. / 6561, -25360. / 2187, 64448. / 6561, -212. / 729},
         {9017. / 3168, -355. / 33, 46732. / 5247, 49. / 176, -5103. / 18656},
         {35. / 384, 0., 500. / 1113, 125. / 192, -2187. / 6784, 11. / 84}}});
};

} // namespace rk
} // namespace Honeycomb

#endif // HC2_RUNGE_KUTTA_HPP