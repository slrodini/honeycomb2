#ifndef DISCRETIZATION_HPP
#define DISCRETIZATION_HPP

#include <honeycomb2/default.hpp>
#include <honeycomb2/utilities.hpp>

#include <honeycomb2/Eigen/Core>

namespace Honeycomb
{

struct StandardGrid {
public:
   StandardGrid(size_t N = 3);

   double interpolate(double t, Eigen::VectorXd &fj, long int start, long int end) const;

   double poli_weight(double t, size_t j) const;
   double poli_weight_sub(double t, size_t j) const;
   double poli_weight_der(double t, size_t j) const;

   double tj(size_t j) const
   {
      return _tj[j];
   }

   size_t _N;
   std::vector<double> _tj;
   std::vector<double> _betaj;
   std::vector<std::vector<double>> _Dij;
};

// ================================================================

// For either radius or angle
struct SingleDiscretizationInfo {
public:
   // Empty default constructor
   SingleDiscretizationInfo()
       : is_periodic(false), to_inter_space(std::function<double(double)>([](double x) {
            return x;
         })),
         to_inter_space_der(std::function<double(double)>([](double x) {
            return x;
         })),
         to_phys_space(std::function<double(double)>([](double x) {
            return x;
         })),
         to_phys_space_der(std::function<double(double)>([](double x) {
            return x;
         })) {};

   SingleDiscretizationInfo(
       std::vector<double> inter, std::vector<size_t> g_size, bool is_per = false,
       std::function<double(double)> to_i_space =
           [](double x) {
              return x;
           },
       std::function<double(double)> to_i_space_der =
           [](double x) {
              (void)x;
              return 1;
           },
       std::function<double(double)> to_p_space =
           [](double x) {
              return x;
           },
       std::function<double(double)> to_p_space_der =
           [](double x) {
              (void)x;
              return 1;
           })
       : is_periodic(is_per), intervals(inter.size() - 1, {0, 0}), intervals_phys(inter.size() - 1, {0, 0}),
         grid_sizes(g_size), to_inter_space(to_i_space), to_inter_space_der(to_i_space_der),
         to_phys_space(to_p_space), to_phys_space_der(to_p_space_der)
   {
      if (g_size.size() != (inter.size() - 1)) {
         logger(Logger::ERROR, std::format("[SingleDiscretizationInfo] Incompatible sizes for number of "
                                           "subintervals ({:d}) andentries in the "
                                           "vector of number of points for each subinterval ({:d})",
                                           inter.size() - 1, g_size.size()));
      }
      for (size_t i = 0; i < inter.size() - 1; i++) {
         intervals[i]      = {to_i_space(inter[i]), to_i_space(inter[i + 1])};
         intervals_phys[i] = {inter[i], inter[i + 1]};
      }
   };

   const bool is_periodic;
   std::vector<std::pair<double, double>> intervals;
   std::vector<std::pair<double, double>> intervals_phys;
   const std::vector<size_t> grid_sizes;
   const std::function<double(double)> to_inter_space;
   const std::function<double(double)> to_inter_space_der;
   const std::function<double(double)> to_phys_space;
   const std::function<double(double)> to_phys_space_der;
};

// Store the grid in either the radius or the angle
struct Grid {

   // u=a <=> t=+1
   // u=b <=> t=-1
   double from_ab_to_m1p1(double u, double a, double b) const noexcept
   {
      return -2 * (u - a) / (b - a) + 1;
   };
   double from_m1p1_to_ab(double t, double a, double b) const noexcept
   {
      return (b - a) * (1 - t) * 0.5 + a;
   };
   double from_ab_to_m1p1_der(double a, double b) const noexcept
   {
      return -2 / (b - a);
   };

   enum W_CASE { N = 1, D = 2, S = 3 };

   Grid(const SingleDiscretizationInfo &d_info);

   // Empty default constructor
   Grid() {};

   double get_der_matrix(size_t a, size_t j, size_t b, size_t k) const;

   // NOTE: u is in interpolation space, not physical space
   template <int w_case>
   double weight_aj(double u, size_t a, size_t j);

   // NOTE: Support is given in interpolation space, not physical space.
   std::pair<double, double> get_support_weight_aj(size_t index) const
   {
      for (size_t a = 0; a < _d_info.intervals.size(); a++) {
         if (_delim_indexes[a] <= index && index < _delim_indexes[a + 1])
            return {_d_info.intervals[a].first, _d_info.intervals[a].second};
      }
      logger(Logger::ERROR, std::format("[Grid::get_support_weight_aj] Index {:d} out of bound", index));
      return {NAN, NAN};
   }

   std::pair<double, double> get_phys_support_weight_aj(size_t index) const
   {
      for (size_t a = 0; a < _d_info.intervals.size(); a++) {
         if (_delim_indexes[a] <= index && index < _delim_indexes[a + 1])
            return {_d_info.intervals_phys[a].first, _d_info.intervals_phys[a].second};
      }
      logger(Logger::ERROR, std::format("[Grid::get_phys_support_weight_aj] Index {:d} out of bound", index));
      return {NAN, NAN};
   }

   const SingleDiscretizationInfo _d_info;
   std::vector<std::function<double(double)>> _weights;
   std::vector<std::function<double(double)>> _weights_der; //  these are dw/du
   std::vector<std::function<double(double)>> _weights_sub; //  these are w-1
   std::vector<double> _coord;
   std::vector<long int> _from_iw_to_ic;
   std::vector<std::vector<long int>> _from_ic_to_iw;
   std::vector<std::vector<double>> _der_matrix;
   std::vector<size_t> _delim_indexes;
   size_t size;
   long int size_li;
   size_t c_size;
   long int c_size_li;
};

namespace RnC
{
struct Pair {
   Pair() : v({0, 0}) {};

   Pair(double a, double b) : v({a, b}) {};

   double &operator[](int i)
   {
      assert(i >= 0 && i <= 1);
      return v[i];
   }
   double operator()(int i) const
   {
      assert(i >= 0 && i <= 1);
      return v[i];
   }
   std::array<double, 2> v;
};

struct Triplet {
   Triplet() : v({0, 0, 0}) {};
   Triplet(double a, double b, double c) : v({a, b, c}) {};

   double &operator[](int i)
   {
      assert(i >= 0 && i <= 2);
      return v[i];
   }
   double operator()(int i) const
   {
      assert(i >= 0 && i <= 2);
      return v[i];
   }

   double distance(const Triplet &other) const
   {
      return sq(v[0] - other.v[0]) + sq(v[1] - other.v[1]) + sq(v[2] - other.v[2]);
   }

   std::array<double, 3> v;

   template <class Archive>
   void serialize(Archive &archive)
   {
      archive(v[0], v[1], v[2]);
   }
};

Triplet from_rhophi_to_x123(double rho, double phi);
Triplet from_rhophi_to_x123(Pair rf);
Pair from_x123_to_rhophi(double x1, double x2, double x3);
Pair from_x123_to_rhophi(Triplet x123);

std::pair<Triplet, Triplet> dx123_drhophi(Pair rhophi);
}; // namespace RnC

struct Grid2D {

   Grid2D(const SingleDiscretizationInfo &d_info_rho, const SingleDiscretizationInfo &d_info_phi);

   // Empty default constructor
   Grid2D()
   {
      is_compliant = false;
   };

   std::function<double(const RnC::Pair &rhophi)> get_dw_dx3_fixed_x1(size_t index) const;

   // NOTE: for the weights
   size_t get_flatten_index(size_t i_r, size_t i_a) const
   {
      return i_r + (grid_radius.size) * i_a;
   }
   size_t get_flatten_index(const std::pair<size_t, size_t> &p) const
   {
      return get_flatten_index(p.first, p.second);
   }
   std::pair<size_t, size_t> get_double_index(size_t index) const
   {
      size_t i_r = index % grid_radius.size;
      size_t i_a = (index - i_r) / grid_radius.size;
      return {i_r, i_a};
   }

   // NOTE: for _x123
   size_t c_get_flatten_index(size_t i_r, size_t i_a) const
   {
      return i_r + (grid_radius.c_size) * i_a;
   }
   size_t c_get_flatten_index(const std::pair<size_t, size_t> &p) const
   {
      return c_get_flatten_index(p.first, p.second);
   }
   std::pair<size_t, size_t> c_get_double_index(size_t index) const
   {
      size_t i_r = index % grid_radius.c_size;
      size_t i_a = (index - i_r) / grid_radius.c_size;
      return {i_r, i_a};
   }

   const Grid grid_radius;
   const Grid grid_angle;

   size_t size;      // Global flatten size
   long int size_li; // Global flatten size

   size_t c_size;      // Size of _fj (different from size, because merging points are stored once)
   long int c_size_li; // Size of _fj (different from size, because merging points are stored once)

   // c_size cooridnates
   std::vector<RnC::Triplet> _x123;

   // Of "size" size, not c_size, it is one for each weight
   std::vector<std::pair<RnC::Triplet, RnC::Triplet>> _x123_minmax;

   // size weights derivatives
   std::vector<std::function<double(const RnC::Pair &rhophi)>> _w;
   std::vector<std::function<double(const RnC::Pair &rhophi)>> _dw_dx3;

   bool is_compliant;

   template <class Archive>
   void serialize(Archive &archive)
   {
      archive(size, c_size);
      archive(_x123);
      archive(_x123_minmax);
   }
};

struct Discretization {
public:
   Discretization(const Grid2D &grid);

   double interpolate_as_weights(const RnC::Pair &rhophi, const Eigen::VectorXd &_fj) const;

   // This only works with sector compliant grids. The sectors are
   //
   //   (++-) -> 0 <= phi <  1
   //   (-+-) -> 1 <= phi <  2
   //   (-++) -> 2 <= phi <  3
   //   (--+) -> 3 <= phi <  4
   //   (+-+) -> 4 <= phi <  5
   //   (+--) -> 5 <= phi <= 6
   //
   double interpolate_as_weights_v2(const RnC::Pair &rhophi, const Eigen::VectorXd &_fj) const;
   double interpolate_as_weights_v3(const RnC::Pair &rhophi, const Eigen::VectorXd &_fj) const;

   double interpolate_as_weights(const RnC::Triplet &x123, const Eigen::VectorXd &_fj) const
   {
      return interpolate_as_weights(RnC::from_x123_to_rhophi(x123), _fj);
   };
   double interpolate_df_dx3_fixed_x1(const RnC::Pair &rhophi, const Eigen::VectorXd &_fj) const;

   Eigen::VectorXd operator()(std::function<double(double, double, double)> const &function) const;

   Eigen::VectorXd discretize(std::function<double(double, double, double)> const &function) const
   {
      return this->operator()(function);
   }

   const Grid2D &_grid;
};

struct LinearGrid {
   static double to_i_space(double x)
   {
      return x;
   };
   static double to_i_space_der(double x)
   {
      (void)x;
      return 1;
   };
   static double to_p_space(double u)
   {
      return u;
   };
   static double to_p_space_der(double u)
   {
      (void)u;
      return 1;
   };
};

struct LogGrid {
   static double to_i_space(double x)
   {
      return log(x);
   };
   static double to_i_space_der(double x)
   {
      return 1.0 / x;
   };
   static double to_p_space(double u)
   {
      return exp(u);
   };
   static double to_p_space_der(double u)
   {
      return exp(u);
   };
};

struct OneOverXGrid {
   static double to_i_space(double x)
   {
      return 1.0 - 1.0 / x;
   };
   static double to_i_space_der(double x)
   {
      return 1.0 / (x * x);
   };
   static double to_p_space(double u)
   {
      return 1.0 / (1.0 - u);
   };
   static double to_p_space_der(double u)
   {
      return 1.0 / sq(1.0 - u);
   };
};

struct OneOverSqrtXGrid {
   static double to_i_space(double x)
   {
      return 1.0 - 1.0 / sqrt(x);
   };
   static double to_i_space_der(double x)
   {
      return 0.5 / sqrt(cu(x));
   };
   static double to_p_space(double u)
   {
      return 1.0 / sq(1.0 - u);
   };
   static double to_p_space_der(double u)
   {
      return 2.0 / cu(1.0 - u);
   };
};

// Utility function, ensures that the angular grid is generated with only 6
// subintervals, corresponding to the 6 triangular regions. This is fundamental for the correct working of
// various parts of the library, especially the routines to determine the support of the weights int the
// physical momentum coordinates. Moreover, the radial grid is defaulted to the logarithmic grid.
Grid2D
generate_compliant_Grid2D(size_t n_pts_for_angle_sector, std::vector<double> radius_inter,
                          std::vector<size_t> radius_g_size,
                          std::function<double(double)> radius_to_i_space     = LogGrid::to_i_space,
                          std::function<double(double)> radius_to_i_space_der = LogGrid::to_i_space_der,
                          std::function<double(double)> radius_to_p_space     = LogGrid::to_p_space,
                          std::function<double(double)> radius_to_p_space_der = LogGrid::to_p_space_der,
                          std::function<double(double)> angle_to_i_space      = LinearGrid::to_i_space,
                          std::function<double(double)> angle_to_i_space_der  = LinearGrid::to_i_space_der,
                          std::function<double(double)> angle_to_p_space      = LinearGrid::to_p_space,
                          std::function<double(double)> angle_to_p_space_der  = LinearGrid::to_p_space_der);

// TODO: Legacy stuff, will be removed later.
struct Discretization1D {
public:
   Discretization1D(const Grid &grid);

   Eigen::VectorXd discretize(const std::function<double(double)> &function) const
   {
      return this->operator()(function);
   }

   double interpolate_as_weights(double x, const Eigen::VectorXd &_fj) const;
   double interpolate_der_as_weights(double x, const Eigen::VectorXd &_fj) const;

   Eigen::VectorXd operator()(const std::function<double(double)> &function) const;

   const Grid &_grid;
};

} // namespace Honeycomb
#endif // DISCRETIZATION_HPP