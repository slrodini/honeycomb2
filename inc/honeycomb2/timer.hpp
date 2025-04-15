#ifndef TIMER_HPP
#define TIMER_HPP

#include <chrono>
#include <string>
#include <iostream>

namespace Honeycomb
{
namespace timer
{
using mark = std::chrono::time_point<std::chrono::high_resolution_clock>;

inline mark now()
{
   return std::chrono::high_resolution_clock::now();
}

double elapsed_ns(const mark &end, const mark &begin)
{
   return static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count());
}

double elapsed_us(const mark &end, const mark &begin)
{
   return 1.0e-3 * static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count());
}

double elapsed_ms(const mark &end, const mark &begin)
{
   return 1.0e-6 * static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count());
}

} // namespace timer
} // namespace Honeycomb

#endif // TIMER_HPP