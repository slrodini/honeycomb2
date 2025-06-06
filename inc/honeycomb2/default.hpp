#ifndef HC2_DEFAULT_HPP
#define HC2_DEFAULT_HPP

/**
 * @file default.hpp
 * @author Simone Rodini (rodini.simone.luigi@gmail.com)
 * @brief  Default includes across Honeycomb
 * @version 0.1
 * @date 2025-06-04
 *
 * @copyright Copyright (c) 2025
 *
 */

// Input/output
#include <cstdio>
#include <iostream>
#include <format>
#include <filesystem>
#include <fstream>

// Memory
#include <cstdlib>

// Math and numerical types
#include <cmath>
#include <cstdint>

// Data structures
#include <map>
#include <vector>
#include <string>
#include <array>
#include <span>
#include <tuple>

// Routines
#include <functional>
#include <algorithm>
#include <random>

// Concepts
#include <type_traits>
#include <concepts>

/// @brief Compure the square of a value
inline double sq(double x) noexcept
{
   return x * x;
};

/// @brief Compure the cube of a value
inline double cu(double x) noexcept
{
   return x * x * x;
};

#endif // HC2_DEFAULT_HPP