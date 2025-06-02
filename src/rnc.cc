#include <honeycomb2/discretization.hpp>

namespace Honeycomb
{
namespace RnC
{
Triplet from_rhophi_to_x123(double rho, double phi)
{
   Triplet x123(0, 0, 0);
   if (phi < 0 || phi > 6) {
      logger(Logger::ERROR,
             std::format("[Grid2D::from_rhophi_to_x123] Invalid angle: {:+.16e}", phi));
   } else if (0 <= phi && phi < 1) {
      x123[0] = rho * (1 - phi);
      x123[1] = rho * (phi);
      x123[2] = rho * (-1);
   } else if (1 <= phi && phi < 2) {
      x123[0] = rho * (1 - phi);
      x123[1] = rho * (1);
      x123[2] = rho * (phi - 2);
   } else if (2 <= phi && phi < 3) {
      x123[0] = rho * (-1);
      x123[1] = rho * (3 - phi);
      x123[2] = rho * (phi - 2);
   } else if (3 <= phi && phi < 4) {
      x123[0] = rho * (phi - 4);
      x123[1] = rho * (3 - phi);
      x123[2] = rho * (1);
   } else if (4 <= phi && phi < 5) {
      x123[0] = rho * (phi - 4);
      x123[1] = rho * (-1);
      x123[2] = rho * (5 - phi);
   } else if (5 <= phi && phi <= 6) {
      x123[0] = rho * (1);
      x123[1] = rho * (phi - 6);
      x123[2] = rho * (5 - phi);
   }
   return x123;
}

Pair from_x123_to_rhophi(double x1, double x2, double x3)
{
   Pair rf(0, 0);
   const double r = std::max(std::fabs(x1), std::max(std::fabs(x2), std::fabs(x3)));

   rf[0] = r;

   if (x1 > 0 && x2 >= 0 && x3 < 0) rf[1] = x2 / r;
   else if (x1 <= 0 && x2 > 0 && x3 < 0) rf[1] = (1 - x1 / r);
   else if (x1 < 0 && x2 > 0 && x3 >= 0) rf[1] = (3 - x2 / r);
   else if (x1 < 0 && x2 <= 0 && x3 > 0) rf[1] = (3 - x2 / r);
   else if (x1 >= 0 && x2 < 0 && x3 > 0) rf[1] = (4 + x1 / r);
   else if (x1 > 0 && x2 < 0 && x3 <= 0) rf[1] = (6 + x2 / r);
   return rf;
}

Triplet from_rhophi_to_x123(Pair rf)
{
   return from_rhophi_to_x123(rf[0], rf[1]);
};
Pair from_x123_to_rhophi(Triplet x123)
{
   return from_x123_to_rhophi(x123[0], x123[1], x123[2]);
};

std::pair<Triplet, Triplet> dx123_drhophi(Pair rhophi)
{
   Triplet dx123_drho, dx123_dphi;

   const double rho = rhophi(0);
   const double phi = rhophi(1);

   if (phi < 0 || phi > 6) {
      logger(Logger::ERROR, std::format("[Grid2D::dx123_drhophi] Invalid angle: {:+.16e}", phi));
   } else if (0 <= phi && phi < 1) {
      dx123_drho[0] = (1 - phi);
      dx123_drho[1] = (phi);
      dx123_drho[2] = (-1);

      dx123_dphi[0] = rho * (-1);
      dx123_dphi[1] = rho * (+1);
      dx123_dphi[2] = rho * (+0);
   } else if (1 <= phi && phi < 2) {
      dx123_drho[0] = (1 - phi);
      dx123_drho[1] = (1);
      dx123_drho[2] = (phi - 2);

      dx123_dphi[0] = rho * (-1);
      dx123_dphi[1] = rho * (+0);
      dx123_dphi[2] = rho * (+1);
   } else if (2 <= phi && phi < 3) {
      dx123_drho[0] = (-1);
      dx123_drho[1] = (3 - phi);
      dx123_drho[2] = (phi - 2);

      dx123_dphi[0] = rho * (+0);
      dx123_dphi[1] = rho * (-1);
      dx123_dphi[2] = rho * (+1);
   } else if (3 <= phi && phi < 4) {
      dx123_drho[0] = (phi - 4);
      dx123_drho[1] = (3 - phi);
      dx123_drho[2] = (1);

      dx123_dphi[0] = rho * (+1);
      dx123_dphi[1] = rho * (-1);
      dx123_dphi[2] = rho * (+0);
   } else if (4 <= phi && phi < 5) {
      dx123_drho[0] = (phi - 4);
      dx123_drho[1] = (-1);
      dx123_drho[2] = (5 - phi);

      dx123_dphi[0] = rho * (+1);
      dx123_dphi[1] = rho * (+0);
      dx123_dphi[2] = rho * (-1);
   } else if (5 <= phi && phi <= 6) {
      dx123_drho[0] = (1);
      dx123_drho[1] = (phi - 6);
      dx123_drho[2] = (5 - phi);

      dx123_dphi[0] = rho * (+0);
      dx123_dphi[1] = rho * (+1);
      dx123_dphi[2] = rho * (-1);
   }

   return {dx123_drho, dx123_dphi};
}
} // namespace RnC
} // namespace Honeycomb