#ifndef HONEYCOMB2_ALPHA_S_HPP
#define HONEYCOMB2_ALPHA_S_HPP

#include <honeycomb2/default.hpp>

/**
 * @file alpha_s.hpp
 * @author Simone Rodini (rodini.simone.luigi@gmail.com)
 * @brief  Get function to compute \f$ a_s(log(Q^2))\f$
 * @version 0.1
 * @date 2025-06-04
 *
 * @copyright Copyright (c) 2025
 *
 */

namespace Honeycomb
{
/**
 * \brief Returns a function to compute \f$ a_s(log(Q^2)) = \alpha_s(log(Q^2)) / 4\pi\f$.
 * \param thresholds   The mass squared of the quarks (in \f$ GeV^2\f$ units).
 * \param Q2ref        The reference scale to fix the value of \f$ \alpha_s \f$ in \f$ GeV^2\f$.
 *                     (default: \f$ M_Z^2 = (91.188)^2 \text{ GeV}^2\f$ ).
 * \param alpha_s_ref  The value of \f$ \alpha_s(Q^2_{ref}) \f$ (default: 0.118)
 * \param Nc           The number of colors. (default: 3)
 */
std::function<double(double)> GetAlphaS_o_4pi(std::array<double, 6> thresholds,
                                              double Q2ref       = 8315.251344,
                                              double alpha_s_ref = 0.118, double Nc = 3);
} // namespace Honeycomb
#endif
