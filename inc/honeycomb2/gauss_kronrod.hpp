#ifndef HC2_GAUSS_KRONROD_HPP
#define HC2_GAUSS_KRONROD_HPP

#include <honeycomb2/default.hpp>

namespace Honeycomb
{
namespace Integration
{

template <typename Rule, std::size_t N>
concept IsGKRule = requires(Rule r) {
   { r.size } -> std::convertible_to<size_t>;
   requires std::same_as<decltype(r.x_gk), const std::array<double, N>>;
   requires std::same_as<decltype(r.w_g), const std::array<double, (N - 1) / 2>>;
   requires std::same_as<decltype(r.w_gk), const std::array<double, N>>;
};

struct GK_21 {
   constexpr static size_t size = 11;
   const static std::array<double, 11> x_gk;
   const static std::array<double, 5> w_g;
   const static std::array<double, 11> w_gk;
};

struct GK_41 {
   constexpr static size_t size = 21;
   const static std::array<double, 21> x_gk;
   const static std::array<double, 10> w_g;
   const static std::array<double, 21> w_gk;
};

struct GK_61 {
   constexpr static size_t size = 31;
   const static std::array<double, 31> x_gk;
   const static std::array<double, 15> w_g;
   const static std::array<double, 31> w_gk;
};

template <typename Rule>
requires IsGKRule<Rule, Rule::size>
struct GaussKronrod {

   // Note: here I do not return the information on the error outside.
   static double integrate(std::function<double(double)> const &fnc, double const &a, double const &b,
                           double eps_rel = 1.0e-10, double eps_abs = 1.0e-10);

private:
   static std::tuple<double, double, double>
   gauss_kronrod_recursive_step(std::function<double(double)> const &fnc, double const &a, double const &b,
                                size_t depth, double eps_rel, double eps_abs);
   static std::tuple<double, double, double>
   gauss_kronrod_simplified(std::function<double(double)> const &fnc, double const &a, double const &b);
};

} // namespace Integration
} // namespace Honeycomb

#endif // HC2_GAUSS_KRONROD_HPP