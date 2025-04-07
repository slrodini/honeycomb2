#pragma once

#include <honeycomb2/discretization.hpp>

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
};

// G G kernels

struct Hhat12GG {
};

struct Hhat31GG {
};

struct Hhat23GG {
};

struct Hplus12GG {
};

struct Hplus13GG {
};

struct Htildeplus12GG {
};

struct Htildeplus13GG {
};

struct Hminus12GG {
};

struct Hminus13GG {
};

struct Vplus13 {
};

struct Vminus13 {
};

struct Wplus13 {
};

struct Wminus13 {
};

struct DeltaW {
};
} // namespace Honeycomb