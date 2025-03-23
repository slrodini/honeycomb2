#include <honeycomb2/honeycomb2.hpp>

int main()
{
   std::function<double(double)> fnc = [](double x) -> double {
      return x * x * exp(x);
   };

   using integrator = Honeycomb::Integration::GaussKronrod<Honeycomb::Integration::GK_41>;
   Honeycomb::logger(Honeycomb::Logger::INFO, std::format("{:.16f}", integrator::integrate(fnc, 0, 1).first));
}