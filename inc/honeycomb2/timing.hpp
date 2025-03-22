#pragma once

#include <chrono>
#include <string>
#include <iostream>

namespace timing
{

using time_point = std::chrono::time_point<std::chrono::high_resolution_clock>;
inline time_point now()
{
   return std::chrono::high_resolution_clock::now();
}

using nanoseconds  = std::chrono::nanoseconds;
using microseconds = std::chrono::microseconds;
using milliseconds = std::chrono::milliseconds;

template <typename T>
struct TimeUnit {
   static std::string get()
   {
      return "Unknown time unit";
   }
};

template <>
struct TimeUnit<nanoseconds> {
   static std::string get()
   {
      return "nanoseconds";
   }
};

template <>
struct TimeUnit<microseconds> {
   static std::string get()
   {
      return "microseconds";
   }
};

template <>
struct TimeUnit<milliseconds> {
   static std::string get()
   {
      return "milliseconds";
   }
};

template <typename time_unit>
double runtime(const time_point &start, const time_point &stop)
{
   return static_cast<double>(std::chrono::duration_cast<time_unit>(stop - start).count());
}

} // namespace timing
