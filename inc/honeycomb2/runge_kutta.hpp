#ifndef RUNGE_KUTTA_HPP
#define RUNGE_KUTTA_HPP

#include <functional>
#include <concepts>
#include <type_traits>
#include <cmath>

#include <iostream>
#include <vector>
#include <array>

namespace Honeycomb
{
namespace runge_kutta
{
// To ensure that Kernel and Solution have all the required operations correctly defined

// Kernel * Solution
template <typename K, typename S>
concept KernelSolutionMultiplicable = requires(K k, S s) {
   { k *s } -> std::convertible_to<S>;
};

// Solution * (any arithmetic type, i.e. anything representable as a pure number)
// and (arithemtic) * Solution
// and Solution / (arithemtic)
template <typename S, typename T>
concept SolutionScalarMultiplicable = std::is_arithmetic_v<T> && requires(S s, T d) {
   { s *d } -> std::convertible_to<S>;
   { s / d } -> std::convertible_to<S>;
   { d *s } -> std::convertible_to<S>;
};

// Solution += Solution
template <typename S>
concept SolutionAddAssign = requires(S s1, S s2) {
   { s1 += s2 } -> std::convertible_to<S &>;
};

// Combined concept, a bit long but list the fundamental types accepted
template <typename K, typename S>
concept RungeKuttaCompatible_1 =
    KernelSolutionMultiplicable<K, S> &&
    (SolutionScalarMultiplicable<S, int> && SolutionScalarMultiplicable<S, float> &&
     SolutionScalarMultiplicable<S, double> && SolutionScalarMultiplicable<S, long double>) &&
    SolutionAddAssign<S>;

template <typename K, typename S>
concept RungeKuttaCompatible_2 = requires {
   // Exact match for: void _copy(const S&)
   { static_cast<void (S::*)(const S &)>(&S::_copy) };

   // Exact match for: void _plus_eq(double, const S&)
   { static_cast<void (S::*)(double, const S &)>(&S::_plus_eq) };

   // Exact match for: void _ker_mul(double, const K&)
   { static_cast<void (S::*)(double, const K &)>(&S::_ker_mul) };
};

// template <typename K, typename S>
// concept RungeKuttaCompatible_2 = requires(S a, const S &b, double x, const K &k) {
//    // This should implement a = b
//    { a._copy(b) } -> std::same_as<void>;

//    // This shloud implement a += x * b
//    // { a._plus_eq(x, b) } -> std::same_as<void>;
//    requires requires {
//       { std::declval<S>()._plus_eq(std::declval<double>(), std::declval<const S &>()) } -> std::same_as<void>;
//    };

//    // This should implement a = x * (k * a)
//    // { a._ker_mul(x, k) } -> std::same_as<void>;
//    requires requires {
//       { std::declval<S>()._ker_mul(std::declval<double>(), std::declval<const K &>()) } -> std::same_as<void>;
//    };
// };

template <typename K, typename S>
concept RungeKuttaCompatible = RungeKuttaCompatible_1<K, S> || RungeKuttaCompatible_2<K, S>;

template <int Order>
struct GenericRungeKuttaTableaux {
   GenericRungeKuttaTableaux() = delete;

   constexpr GenericRungeKuttaTableaux(std::array<double, Order> const &ci, std::array<double, Order> const &bi,
                                       std::array<std::array<double, Order>, Order> const &ai)
       : _ai(ai), _bi(bi), _ci(ci) {};

   const std::array<std::array<double, Order>, Order> _ai;
   const std::array<double, Order> _bi;
   const std::array<double, Order> _ci;
};

template <typename Kernel, typename Solution, int Order>
requires RungeKuttaCompatible<Kernel, Solution>
class GenericRungeKutta
{
public:
   // Templetized over arguments for Solution constructor
   template <typename... Args>
   GenericRungeKutta(
       Kernel const &K, Solution const &S, GenericRungeKuttaTableaux<Order> const &tab,
       std::function<double(double)> as =
           [](double t) {
              (void)t;
              return 1;
           },
       double pref = 1, double t = 0, double dt = 0.01,
       std::function<void(double, Kernel &, Solution &)> c =
           [](double t, Kernel &K, Solution &S) {
              (void)t;
              (void)K;
              (void)S;
              return;
           },
       Args &&...args)
       : _t(t), _dt(dt), _prefactor(pref), _callback(std::move(c)), _As_t(std::move(as)), _solution(S), _temp_step(S),
         _kernel(K), _ai(tab._ai), _bi(tab._bi), _ci(tab._ci), callback_each_step(true)
   {
      for (auto &ki : _ki)
         ki = Solution(std::forward<Args>(args)...);
   };

   // Suggestion: overload operator= and operator += in solution to avoid memory re-allocation every time
   void operator()()
   {
      if constexpr (RungeKuttaCompatible_1<Kernel, Solution>) {
         _temp_step = _solution;

         for (int i = 0; i < Order; i++) {

            _ki[i] = _solution;

            for (int j = 0; j < i; j++)
               _ki[i] += _ai[i][j] * _ki[j];

            const double pref = _prefactor * _dt * _As_t(_t + _dt * _ci[i]);
            _ki[i]            = pref * (_kernel * _ki[i]);

            _temp_step += _bi[i] * _ki[i];
         }

         _solution = _temp_step;
      } else {

         _temp_step._copy(_solution);

         for (int i = 0; i < Order; i++) {

            _ki[i]._copy(_solution);

            for (int j = 0; j < i; j++)
               _ki[i]._plus_eq(_ai[i][j], _ki[j]);

            const double pref = _prefactor * _dt * _As_t(_t + _dt * _ci[i]);

            _ki[i]._ker_mul(pref, _kernel);

            _temp_step._plus_eq(_bi[i], _ki[i]);
         }

         _solution._copy(_temp_step);
      }

      _t += _dt;
      if (callback_each_step) _callback(_t, _kernel, _solution);
   };

   void operator()(std::vector<double> const &thresholds, size_t n_step = 10)
   {
      // Remove callback at each step
      callback_each_step = false;
      for (const double &tf : thresholds) {
         _dt = (tf - _t) / static_cast<double>(n_step);
         for (size_t i = 0; i < n_step; i++)
            (*this)();
         if (std::fabs(_t - tf) > 1.0e-14) {
            std::cerr << "RungeKutta Error: not reached the correct scale. This is a bug. Aborting\n";
            exit(-1);
         }
         _t = tf; // Set exactly equal to threshold
         _callback(tf, _kernel, _solution);
      }
   }

   void reset(Solution const &S, double t)
   {
      _t                 = t;
      callback_each_step = true;
      _solution          = S;
   }

   Solution GetSolution() const
   {
      return _solution;
   }
   Kernel GetKernel() const
   {
      return _kernel;
   }

   void SetCallbackFunction(std::function<void(double, Kernel &, Solution &)> c)
   {
      _callback = c;
   };

private:
   double _t, _dt;    // time, time-step and time-step / 2
   double _prefactor; // eventual prefactor to as(t) * K * S (e.g. 2 or -1 or color factor, depending on definitions)
   std::array<Solution, Order> _ki; // intermediate variables
   std::function<void(double, Kernel &, Solution &)>
       _callback;                       // Callback function for each step, initlialzed to the 'do-nothing' function
   std::function<double(double)> _As_t; // Function to compute as(t) = \alpha_s / 4\pi at given t = \log(\mu^2)

   Solution _solution;
   Solution _temp_step;
   Kernel _kernel;

   std::array<std::array<double, Order>, Order> _ai;
   std::array<double, Order> _bi;
   std::array<double, Order> _ci;
   bool callback_each_step;
};

constexpr GenericRungeKuttaTableaux<4> RKOriginal({0, 0.5, 0.5, 1}, {1.0 / 6.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 6.0},
                                                  {{{0, 0, 0, 0}, {0.5, 0, 0, 0}, {0, 0.5, 0, 0}, {0, 0, 1, 0}}});

constexpr GenericRungeKuttaTableaux<13> DOPRI8(
    {0., 1. / 18, 1. / 12, 1. / 8, 5. / 16, 3. / 8, 59. / 400, 93. / 200, 5490023248. / 9719169821, 13. / 20,
     1201146811. / 1299019798, 1., 1.},
    {14005451. / 335480064, 0., 0., 0., 0., -59238493. / 1068277825, 181606767. / 758867731, 561292985. / 797845732,
     -1041891430. / 1371343529, 760417239. / 1151165299, 118820643. / 751138087, -528747749. / 2220607170, 1. / 4},
    {{{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {1. / 18, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {1. / 48, 1. / 16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {1. / 32, 0., 3. / 32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {
          5. / 16,
          0.,
          -75. / 64,
          75. / 64,
          0,
          0,
          0,
          0,
          0,
          0,
          0,
          0,
          0,
      },
      {
          3. / 80,
          0.,
          0.,
          3. / 16,
          3. / 20,
          0,
          0,
          0,
          0,
          0,
          0,
          0,
          0,
      },
      {29443841. / 614563906, 0., 0., 77736538. / 692538347, -28693883. / 1125000000, 23124283. / 1800000000, 0, 0, 0,
       0, 0, 0, 0},
      {16016141. / 946692911, 0., 0., 61564180. / 158732637, 22789713. / 633445777, 545815736. / 2771057229,
       -180193667. / 1043307555, 0, 0, 0, 0, 0, 0},
      {39632708. / 573591083, 0., 0., -433636366. / 683701615, -421739975. / 2616292301, 100302831. / 723423059,
       790204164. / 839813087, 800635310. / 3783071287, 0, 0, 0, 0, 0},
      {246121993. / 1340847787, 0., 0., -37695042795. / 15268766246, -309121744. / 1061227803, -12992083. / 490766935,
       6005943493. / 2108947869, 393006217. / 1396673457, 123872331. / 1001029789, 0, 0, 0, 0},
      {-1028468189. / 846180014, 0., 0., 8478235783. / 508512852, 1311729495. / 1432422823, -10304129995. / 1701304382,
       -48777925059. / 3047939560, 15336726248. / 1032824649, -45442868181. / 3398467696, 3065993473. / 597172653, 0, 0,
       0},
      {185892177. / 718116043, 0., 0., -3185094517. / 667107341, -477755414. / 1098053517, -703635378. / 230739211,
       5731566787. / 1027545527, 5232866602. / 850066563, -4093664535. / 808688257, 3962137247. / 1805957418,
       65686358. / 487910083, 0, 0},
      {403863854. / 491063109, 0., 0., -5068492393. / 434740067, -411421997. / 543043805, 652783627. / 914296604,
       11173962825. / 925320556, -13158990841. / 6184727034, 3936647629. / 1978049680, -160528059. / 685178525,
       248638103. / 1413531060, 0., 0}}}

);

} // namespace runge_kutta
} // namespace Honeycomb

#endif // RUNGE_KUTTA_HPP