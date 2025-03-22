#include <honeycomb2/random_engine.hpp>

namespace Honeycomb
{
namespace Random
{
std::shared_ptr<RandomEngine> get_random_engine(uint32_t seed)
{
   static std::shared_ptr<Honeycomb::Random::RandomEngine> engine =
       std::make_shared<Honeycomb::Random::RandomEngine>(seed);
   return engine;
}

double random_uniform()
{
   std::uniform_real_distribution<double> dist(0.0, 1.0);
   return dist(*get_random_engine());
}

double random_uniform(double low, double high)
{
   return (high - low) * random_uniform() + low;
}

double random_normal()
{
   std::normal_distribution<double> dist(0.0, 1.0);
   return dist(*get_random_engine());
}

std::vector<double> random_normal(size_t n)
{
   std::vector<double> res(n, 0);
   for (double &v : res)
      v = random_normal();
   return res;
}
} // namespace Random
} // namespace Honeycomb