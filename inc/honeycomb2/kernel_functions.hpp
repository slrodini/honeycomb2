#ifndef KERNEL_FUNCTIONS_HPP
#define KERNEL_FUNCTIONS_HPP

#include <honeycomb2/discretization.hpp>
#include <honeycomb2/utilities.hpp>

namespace Honeycomb
{

// Qbar Q kernels
struct Hhat12 {
   // c_a is _x123 index, in [0, g.c_size)
   // aP is _w index, in [0, g.size)
   static double subtracted_integrate(size_t c_a, size_t aP, const Grid2D &g);
};

struct Hhat13 {
   static double subtracted_integrate(size_t c_a, size_t aP, const Grid2D &g);
};

struct Hhat23 {
   static double subtracted_integrate(size_t c_a, size_t aP, const Grid2D &g);
};

struct Hplus12 {
   static double integrate(size_t c_a, size_t aP, const Grid2D &g);
   // static double integrate_v2(size_t c_a, size_t c_aP, const Grid2D &g);
};

struct Hplus13 {
   static double integrate(size_t c_a, size_t aP, const Grid2D &g);
};

struct Hplus23 {
   static double integrate(size_t c_a, size_t aP, const Grid2D &g);
};

struct He23P23 {
   static double integrate(size_t c_a, size_t aP, const Grid2D &g);
};

struct Hminus12 {
   static double integrate(size_t c_a, size_t aP, const Grid2D &g);
};

struct Hminus23 {
   static double integrate(size_t c_a, size_t aP, const Grid2D &g);
};

struct Hd13 {
   static double integrate(size_t c_a, size_t aP, const Grid2D &g);
};

// G G kernels

struct Hhat12GG {
   static double subtracted_integrate(size_t c_a, size_t aP, const Grid2D &g);
};

struct Hhat23GG {
   static double subtracted_integrate(size_t c_a, size_t aP, const Grid2D &g);
};

struct Hhat31GG {
   static double subtracted_integrate(size_t c_a, size_t aP, const Grid2D &g);
};

// This is already sum of -4 H^+_{12} - 2 \tilde{H}^+_{12}
struct HplusSum12GG {
   static double integrate(size_t c_a, size_t aP, const Grid2D &g);
};

// This is already sum of -4 H^+_{13} - 2 \tilde{H}^+_{13}
struct HplusSum13GG {
   static double integrate(size_t c_a, size_t aP, const Grid2D &g);
};

// These have alreay the factor 6 inside!
struct Hminus12GG {
   static double integrate(size_t c_a, size_t aP, const Grid2D &g);
};

// These have alreay the factor 6 inside!
struct Hminus13GG {
   static double integrate(size_t c_a, size_t aP, const Grid2D &g);
};

// QG and GQ mixing

struct Vplus13 {
   static double integrate(size_t c_a, size_t aP, const Grid2D &g);
   static double subtracted_integrate(size_t c_a, size_t aP, const Grid2D &g);
};

struct Vminus13 {
   static double integrate(size_t c_a, size_t aP, const Grid2D &g);
};

struct Wplus13 {
   static double integrate(size_t c_a, size_t aP, const Grid2D &g);
};

struct Wminus13 {
   static double integrate(size_t c_a, size_t aP, const Grid2D &g);
};

struct DeltaW {
   static double integrate(size_t c_a, size_t aP, const Grid2D &g);
};

struct Wplus13P23 {
   static double integrate(size_t c_a, size_t aP, const Grid2D &g);
};

struct Wminus13P23 {
   static double integrate(size_t c_a, size_t aP, const Grid2D &g);
};

struct DeltaWP23 {
   static double integrate(size_t c_a, size_t aP, const Grid2D &g);
};

} // namespace Honeycomb

#endif // KERNEL_FUNCTIONS_HPP