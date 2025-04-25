#ifndef RANDOM_ENGINE_HPP
#define RANDOM_ENGINE_HPP

#include <honeycomb2/default.hpp>

namespace Honeycomb
{
namespace Random
{
using RandomEngine = std::ranlux48;

std::shared_ptr<RandomEngine> get_random_engine(uint32_t seed = 12345);

// If a seeded engine is need this function must be called AT THE BEGINNING of the main.
inline void seed_random_engine(uint32_t seed)
{
   (void)get_random_engine(seed);
}

double random_uniform();
double random_uniform(double low, double high);
std::vector<double> random_uniform(size_t n, double low, double high);
double random_normal();
std::vector<double> random_normal(size_t n, double mu = 0, double sigma = 1);

} // namespace Random
} // namespace Honeycomb

#endif // RANDOM_ENGINE_HPP