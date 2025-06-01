#ifndef HONEYCOMB2_ALPHA_S_HPP
#define HONEYCOMB2_ALPHA_S_HPP

#include <honeycomb2/default.hpp>
namespace Honeycomb
{
/**
 * \brief Returns a function to compute as(Log(Q^2))
 */
std::function<double(double)> GetAlphaS_o_4pi(std::array<double, 6> thresholds, double Q2ref = 8315.251344,
                                              double alpha_s_ref = 0.118, double Nc = 3);
} // namespace Honeycomb
#endif
