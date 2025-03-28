#pragma once

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

// Routines
#include <functional>
#include <algorithm>
#include <random>

// Concepts
#include <type_traits>
#include <concepts>

// Small inline functions
inline double sq(double x) noexcept
{
   return x * x;
};

inline double cu(double x) noexcept
{
   return x * x * x;
};