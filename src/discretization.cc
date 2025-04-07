#include <honeycomb2/discretization.hpp>

namespace
{
std::map<size_t, Honeycomb::StandardGrid> stored_grids;
}

namespace Honeycomb
{

// =======================================================
// =================== STANDARD GRID =====================
// =======================================================
double linear_map(double x, double a, double b)
{
   return (x - a) / (b - a);
}
double linear_map_inverse(double x, double a, double b)
{
   return (b - a) * x + a;
}

double chebyshev_points_scaled(size_t j, size_t N)
{
   if (j > N) logger(Logger::ERROR, "[Chebyshev_points_scaled]: invalid index.");
   double jd = static_cast<double>(j);
   double Nd = static_cast<double>(N);
   return cos(M_PI * jd / Nd);
}

StandardGrid::StandardGrid(size_t N)
    : _N(N), _tj(N + 1, 0), _betaj(N + 1, 1), _Dij(N + 1, std::vector<double>(N + 1, 0))
{
   _betaj[0] = 0.5;
   _betaj[N] = 0.5;
   for (size_t j = 0; j <= N; j++) {
      _tj[j]     = chebyshev_points_scaled(j, N);
      _betaj[j] *= (j % 2 ? 1 : -1);
   }

   _Dij[0][0] = (2 * N * N + 1) / 6.0;
   _Dij[N][N] = -((2.0 * N * N + 1) / 6.0);
   for (size_t j = 1; j < N; j++) {
      _Dij[j][j] = (-_tj[j]) / (2.0 * (1.0 - _tj[j] * _tj[j]));
   }
   for (size_t j = 0; j <= N; j++) {
      for (size_t k = 0; k < j; k++) {
         _Dij[j][k] = (_betaj[k] / _betaj[j]) / (_tj[j] - _tj[k]);
      }
      for (size_t k = j + 1; k <= N; k++) {
         _Dij[j][k] = (_betaj[k] / _betaj[j]) / (_tj[j] - _tj[k]);
      }
   }
}

double StandardGrid::interpolate(double t, Eigen::VectorXd &fj, long int start, long int end) const
{
   if (t < -1 || t > 1 || (end - start) != static_cast<long int>(_N)) {
      logger(
          Logger::ERROR,
          std::format(
              "[StandardGrid::interpolate]: t={:+.16e} \\notin [-1, +1] OR view into fj of wrong size: [{:+d}, {:+d}]",
              t, start, end));
   }
   double sum = 0;
   size_t j   = 0;
   for (long int i = start; i <= end; i++, j++) {
      sum += poli_weight(t, j) * fj[i];
   }
   return sum;
}

double StandardGrid::poli_weight(double t, size_t j) const
{
   if (t < -1 || t > 1) {
      logger(Logger::ERROR, std::format("[StandardGrid::poli_weight]: t={:+.16e} \\notin [-1, +1]", t));
   }
   if (fabs(t - _tj[j]) < 1.0e-15) return 1.0;

   // Baricentric form
   double den = 0;
   for (size_t l = 0; l <= _N; l++) {
      if (fabs(t - _tj[l]) < 1.0e-15) return 0.0;

      den += _betaj[l] / (t - _tj[l]);
   }
   return _betaj[j] / den / (t - _tj[j]);
}

// bj(u) - 1
double StandardGrid::poli_weight_sub(double t, size_t j) const
{
   if (t < -1 || t > 1) {
      logger(Logger::ERROR, std::format("[StandardGrid::poli_weight_sub]: t={:+.16e} \\notin [-1, +1]", t));
   }
   // Baricentric form
   double den = 0;
   double num = 0;
   for (size_t l = 0; l <= _N; l++) {
      if (fabs(t - _tj[l]) < 1.0e-15) return (l == j) ? 0 : -1;
      den += _betaj[l] / (t - _tj[l]);
      if (l != j) num += _betaj[l] / (t - _tj[l]);
   }
   return -num / den;
}

double StandardGrid::poli_weight_der(double t, size_t j) const
{
   if (t < -1 || t > 1) {
      logger(Logger::ERROR, std::format("[StandardGrid::poli_weight_der]: t={:+.16e} \\notin [-1, +1]", t));
   }
   double res = 0;
   for (size_t i = 0; i <= _N; i++) {
      res += poli_weight(t, i) * _Dij[i][j];
   }
   return res;
}

// =======================================================
// ========================= GRID ========================
// =======================================================

Grid::Grid(const SingleDiscretizationInfo &d_info) : _d_info(d_info)
{
   size = 0;
   for (size_t i = 0; i < d_info.intervals.size(); i++) {
      // NOTE: This works ok because grid_sizes stores N, but the points are N+1 always
      // NOTE: This is implementation of non-overlapping grids, leading to simpler weights
      // TODO: Interpolation was good in tests, but I need to check that the Kernels are not screwed by this
      size += d_info.grid_sizes[i] + 1;
      if (stored_grids.find(d_info.grid_sizes[i]) == stored_grids.end()) {
         stored_grids.emplace(d_info.grid_sizes[i], StandardGrid(d_info.grid_sizes[i]));
      }
   }

   size_li   = static_cast<long int>(size);
   c_size    = size - d_info.intervals.size() + (!d_info.is_periodic);
   c_size_li = static_cast<long int>(size - d_info.intervals.size() + (!d_info.is_periodic));

   _weights.resize(size);
   _weights_der.resize(size);
   _weights_sub.resize(size);
   _from_iw_to_ic.resize(size);
   _coord.resize(c_size);
   _delim_indexes.resize(d_info.intervals.size() + 1);

   _der_matrix.resize(size);
   for (size_t i = 0; i < size; i++)
      _der_matrix[i].resize(size);

   // No check on order of the interval
   long int index       = 0;
   long int index_coord = 0;
   for (size_t a = 0; a < d_info.intervals.size(); a++) {

      const StandardGrid &sg = stored_grids.at(d_info.grid_sizes[a]);
      _delim_indexes[a]      = index;

      for (size_t j = 0; j <= d_info.grid_sizes[a]; j++) {
         // NOTE: _coord stores only the UNIQUES coordinates

         if (j != d_info.grid_sizes[a] || (!d_info.is_periodic && a == d_info.intervals.size() - 1))
            _coord[index_coord] =
                d_info.to_phys_space(from_m1p1_to_ab(sg.tj(j), d_info.intervals[a].first, d_info.intervals[a].second));

         _from_iw_to_ic[index] = index_coord;
         if (j != d_info.grid_sizes[a]) index_coord++;

         _weights[index] = [a, j, this](double u) -> double {
            return this->weight_aj<W_CASE::N>(u, a, j);
         };
         _weights_der[index] = [a, j, this](double u) -> double {
            return this->weight_aj<W_CASE::D>(u, a, j);
         };
         _weights_sub[index] = [a, j, this](double u) -> double {
            return this->weight_aj<W_CASE::S>(u, a, j);
         };

         long int inner_index = 0;
         for (size_t b = 0; b < d_info.intervals.size(); b++) {
            for (size_t k = 0; k <= d_info.grid_sizes[b]; k++) {
               _der_matrix[index][inner_index++] = get_der_matrix(a, j, b, k);
            }
         }
         index++;
      }
   }
   _delim_indexes[d_info.intervals.size()] = index;
   if (d_info.is_periodic) _from_iw_to_ic[size - 1] = 0;

   _from_ic_to_iw.resize(c_size);
   for (size_t i = 0; i < _from_iw_to_ic.size(); i++) {
      _from_ic_to_iw[_from_iw_to_ic[i]].push_back(i);
   }

   // for (const auto &v : _from_ic_to_iw) {
   //    logger(Logger::WARNING, std::format("{:d}", v.size()));
   // }
   // logger(Logger::WARNING, "========================");
}

double Grid::get_der_matrix(size_t a, size_t j, size_t b, size_t k) const
{
   if (a != b) return 0.0;
   StandardGrid &sg = stored_grids.at(_d_info.grid_sizes[a]);
   return sg._Dij[j][k];
}

template <int w_case>
double Grid::weight_aj(double u, size_t a, size_t j)
{
   StandardGrid &sg = stored_grids.at(_d_info.grid_sizes[a]);
   double res       = 0;

   bool condition_x =
       a == _d_info.intervals.size() - 1 ? u <= _d_info.intervals[a].second : u < _d_info.intervals[a].second;

   if (u >= _d_info.intervals[a].first && condition_x) {
      if constexpr (w_case == D) {
         const double dl_dx = from_ab_to_m1p1_der(_d_info.intervals[a].first, _d_info.intervals[a].second);
         const double dw_dl =
             sg.poli_weight_der(from_ab_to_m1p1(u, _d_info.intervals[a].first, _d_info.intervals[a].second), j);
         res += dw_dl * dl_dx;
      } else if constexpr (w_case == N) {
         res += sg.poli_weight(from_ab_to_m1p1(u, _d_info.intervals[a].first, _d_info.intervals[a].second), j);
      } else {
         res += sg.poli_weight_sub(from_ab_to_m1p1(u, _d_info.intervals[a].first, _d_info.intervals[a].second), j);
      }
   }

   return res;
}

// =======================================================
// ======================== GRID2D =======================
// =======================================================

Grid2D::Grid2D(const SingleDiscretizationInfo &d_info_rho, const SingleDiscretizationInfo &d_info_phi)
    : grid_radius(d_info_rho), grid_angle(d_info_phi), size(grid_radius.size * grid_angle.size),
      size_li(static_cast<long int>(size)), c_size(grid_radius.c_size * grid_angle.c_size),
      c_size_li(static_cast<long int>(c_size)), _x123(c_size), _x123_minmax(size), _w(size), _dw_dx3(size)
{
   size_t k = 0;
   size_t j = 0;
   // These are loops over the double w indexes
   for (size_t a = 0; a < grid_angle._d_info.intervals.size(); a++) {
      for (size_t ka = 0; ka <= grid_angle._d_info.grid_sizes[a]; ka++) {

         j = 0;
         for (size_t b = 0; b < grid_radius._d_info.intervals.size(); b++) {

            // NOTE: the support only needs the sub-interval
            auto [phi_min, phi_max] = grid_angle._d_info.intervals_phys[a];
            auto [rho_min, rho_max] = grid_radius._d_info.intervals_phys[b];

            RnC::Triplet pts[4] = {
                RnC::from_rhophi_to_x123(rho_min, phi_min),
                RnC::from_rhophi_to_x123(rho_min, phi_max),
                RnC::from_rhophi_to_x123(rho_max, phi_min),
                RnC::from_rhophi_to_x123(rho_max, phi_max),
            };

            RnC::Triplet xmin(pts[0]), xmax(pts[0]);
            for (size_t q = 1; q < 4; q++) {
               const RnC::Triplet &pt = pts[q];
               for (size_t i_x = 0; i_x < 3; i_x++) {
                  xmin[i_x] = std::min(pt(i_x), xmin(i_x));
                  xmax[i_x] = std::max(pt(i_x), xmax(i_x));
               }
            }

            for (size_t kb = 0; kb <= grid_radius._d_info.grid_sizes[b]; kb++) {
               { // Stuff
                  size_t c_k = grid_angle._from_iw_to_ic[k];
                  size_t c_j = grid_radius._from_iw_to_ic[j];

                  size_t index      = get_flatten_index(j, k);
                  long int index_li = static_cast<long int>(index);

                  size_t f_index      = c_get_flatten_index(c_j, c_k);
                  long int f_index_li = static_cast<long int>(f_index);

                  RnC::Pair rf(grid_radius._coord[c_j], grid_angle._coord[c_k]);
                  _x123[f_index] = RnC::from_rhophi_to_x123(rf);

                  _w[index] = [j, k, this](const RnC::Pair &rhophi) -> double {
                     const double r = grid_radius._d_info.to_inter_space(rhophi(0));
                     const double f = grid_angle._d_info.to_inter_space(rhophi(1));
                     return grid_radius._weights[j](r) * grid_angle._weights[k](f);
                  };

                  _dw_dx3[index] = get_dw_dx3_fixed_x1(index);

                  _x123_minmax[index] = {xmin, xmax};
               }
               j++;
            }
         }
         k++;
      }
   }
};

std::function<double(const RnC::Pair &rhophi)> Grid2D::get_dw_dx3_fixed_x1(size_t index) const
{
   auto [j, k] = get_double_index(index);

   std::function<double(const RnC::Pair &rhophi)> res = [this, j, k](const RnC::Pair &rhophi) -> double {
      const Grid &radius = grid_radius;
      const Grid &angle  = grid_angle;

      const double r = radius._d_info.to_inter_space(rhophi(0));
      const double f = angle._d_info.to_inter_space(rhophi(1));

      auto f_supp = angle.get_support_weight_aj(k);
      if (f < f_supp.first || f > f_supp.second) return 0;
      auto r_supp = radius.get_support_weight_aj(j);
      if (r < r_supp.first || r > r_supp.second) return 0;

      const double drho_dr     = radius._d_info.to_phys_space_der(r);
      const double dphi_df     = angle._d_info.to_phys_space_der(f);
      const double drrho_dfphi = drho_dr * dphi_df;

      const auto [dxi_drho, dxi_dphi] = RnC::dx123_drhophi(rhophi);

      const double den = drrho_dfphi * (-dxi_dphi(2) * dxi_drho(0) + dxi_dphi(0) * dxi_drho(2));

      if (std::fabs(den) < 1.0e-15) {
         logger(Logger::ERROR, "[interpolate_df_dx3_fixed_x1] Vanishing denominator, can this be possible?");
      }

      const double num = (angle._weights[k](f) * radius._weights_der[j](r) * (dphi_df * dxi_dphi(0)) -
                          radius._weights[j](r) * angle._weights_der[k](f) * (drho_dr * dxi_drho(0)));

      return num / den;
   };
   return res;
}

// =======================================================
// =================== DISCRETIZATION ====================
// =======================================================

Discretization::Discretization(const Grid2D &grid) : _grid(grid)
{
}

Eigen::VectorXd Discretization::operator()(std::function<double(double, double, double)> const &function) const
{
   Eigen::VectorXd _fj(_grid.c_size_li);
   for (long int i = 0; i < _grid.c_size_li; i++) {
      _fj(i) = function(_grid._x123[i](0), _grid._x123[i](1), _grid._x123[i](2));
   }
   return _fj;
}

double Discretization::interpolate_as_weights(const RnC::Pair &rhophi, const Eigen::VectorXd &_fj) const
{
   if (_fj.size() != _grid.c_size_li) {
      logger(
          Logger::ERROR,
          std::format("[interpolate_df_dx3_fixed_x1] size of _fj ({:d}) does not match size stored in the grid ({:d}).",
                      _fj.size(), _grid.c_size_li));
   }
   const Grid &radius = _grid.grid_radius;
   const Grid &angle  = _grid.grid_angle;
   const double r     = radius._d_info.to_inter_space(rhophi(0));
   const double f     = angle._d_info.to_inter_space(rhophi(1));

   double res = 0;

   for (size_t k = 0; k < angle.size; k++) {
      auto f_supp = angle.get_support_weight_aj(k);
      if (f < f_supp.first || f > f_supp.second) continue;

      for (size_t j = 0; j < radius.size; j++) {
         auto r_supp = radius.get_support_weight_aj(j);
         if (r < r_supp.first || r > r_supp.second) continue;

         size_t c_k = angle._from_iw_to_ic[k];
         size_t c_j = radius._from_iw_to_ic[j];

         const double fj  = _fj(static_cast<long int>(_grid.c_get_flatten_index(c_j, c_k)));
         res             += fj * radius._weights[j](r) * angle._weights[k](f);
      }
   }

   return res;
}

double Discretization::interpolate_as_weights_v3(const RnC::Pair &rhophi, const Eigen::VectorXd &_fj) const
{
   if (_fj.size() != _grid.c_size_li) {
      logger(
          Logger::ERROR,
          std::format("[interpolate_as_weights_v3] size of _fj ({:d}) does not match size stored in the grid ({:d}).",
                      _fj.size(), _grid.c_size_li));
   }
   const Grid &radius = _grid.grid_radius;
   const Grid &angle  = _grid.grid_angle;
   const double rho   = rhophi(0);
   const double phi   = rhophi(1);

   double res = 0;

   for (size_t k = 0; k < angle.size; k++) {
      auto f_supp = angle.get_phys_support_weight_aj(k);
      if (phi < f_supp.first || phi > f_supp.second) continue;

      for (size_t j = 0; j < radius.size; j++) {
         auto r_supp = radius.get_phys_support_weight_aj(j);
         if (rho < r_supp.first || rho > r_supp.second) continue;

         size_t c_k = angle._from_iw_to_ic[k];
         size_t c_j = radius._from_iw_to_ic[j];

         size_t index = _grid.get_flatten_index(j, k);

         const double fj  = _fj(static_cast<long int>(_grid.c_get_flatten_index(c_j, c_k)));
         res             += fj * _grid._w[index](rhophi);
      }
   }

   return res;
}

double Discretization::interpolate_as_weights_v2(const RnC::Pair &rhophi, const Eigen::VectorXd &_fj) const
{
   if (_fj.size() != _grid.c_size_li) {
      logger(
          Logger::ERROR,
          std::format("[interpolate_df_dx3_fixed_x1] size of _fj ({:d}) does not match size stored in the grid ({:d}).",
                      _fj.size(), _grid.c_size_li));
   }
   const Grid &radius = _grid.grid_radius;
   const Grid &angle  = _grid.grid_angle;
   const double r     = radius._d_info.to_inter_space(rhophi(0));
   const double f     = angle._d_info.to_inter_space(rhophi(1));

   size_t sector = static_cast<size_t>(floor(rhophi(1)));
   if (sector == 6) sector--;

   // Works only with Grid2D compliant with the angular sector decomposition
   size_t kmin = angle._delim_indexes[sector];
   size_t kmax = angle._delim_indexes[sector + 1];
   size_t jmin = 0;
   size_t jmax = 0;

   for (size_t a = 0; a < radius._d_info.intervals.size(); a++) {
      size_t low  = radius._from_iw_to_ic[radius._delim_indexes[a]];
      size_t high = radius._from_iw_to_ic[radius._delim_indexes[a + 1] - 1];
      if (radius._coord[low] <= rhophi(0) && rhophi(0) < radius._coord[high]) {
         jmin = radius._delim_indexes[a];
         jmax = radius._delim_indexes[a + 1];
         break;
      }
   }

   double res = 0;

   for (size_t k = kmin; k < kmax; k++) {
      for (size_t j = jmin; j < jmax; j++) {
         size_t c_k       = angle._from_iw_to_ic[k];
         size_t c_j       = radius._from_iw_to_ic[j];
         const double fj  = _fj(static_cast<long int>(_grid.c_get_flatten_index(c_j, c_k)));
         res             += fj * radius._weights[j](r) * angle._weights[k](f);
      }
   }

   return res;
}

double Discretization::interpolate_df_dx3_fixed_x1(const RnC::Pair &rhophi, const Eigen::VectorXd &_fj) const
{

   if (_fj.size() != _grid.c_size_li) {
      logger(
          Logger::ERROR,
          std::format("[interpolate_df_dx3_fixed_x1] size of _fj ({:d}) does not match size stored in the grid ({:d}).",
                      _fj.size(), _grid.c_size_li));
   }
   const Grid &radius = _grid.grid_radius;
   const Grid &angle  = _grid.grid_angle;
   const double r     = radius._d_info.to_inter_space(rhophi(0));
   const double f     = angle._d_info.to_inter_space(rhophi(1));

   double res = 0;

   for (size_t k = 0; k < angle.size; k++) {
      auto f_supp = angle.get_support_weight_aj(k);
      if (f < f_supp.first || f > f_supp.second) continue;

      for (size_t j = 0; j < radius.size; j++) {
         auto r_supp = radius.get_support_weight_aj(j);
         if (r < r_supp.first || r > r_supp.second) continue;

         size_t c_k = angle._from_iw_to_ic[k];
         size_t c_j = radius._from_iw_to_ic[j];

         size_t index = _grid.get_flatten_index(j, k);

         const double fj  = _fj(static_cast<long int>(_grid.c_get_flatten_index(c_j, c_k)));
         res             += fj * _grid._dw_dx3[index](rhophi);
      }
   }
   return res;
}

Grid2D generate_compliant_Grid2D(
    size_t n_pts_for_angle_sector, std::vector<double> radius_inter, std::vector<size_t> radius_g_size,
    std::function<double(double)> radius_to_i_space, std::function<double(double)> radius_to_i_space_der,
    std::function<double(double)> radius_to_p_space, std::function<double(double)> radius_to_p_space_der,
    std::function<double(double)> angle_to_i_space, std::function<double(double)> angle_to_i_space_der,
    std::function<double(double)> angle_to_p_space, std::function<double(double)> angle_to_p_space_der)
{

   // SECTION: Sanity checks on the radial intervals
   if (radius_inter.size() == 0) {
      logger(Logger::ERROR, "[generate_complaiant_Grid2D] Empty radial interval.");
   } else if (radius_inter.size() == 1) {
      if (is_near(radius_inter[0], 1)) {
         logger(Logger::ERROR, "[generate_complaiant_Grid2D] Radial interval with only one point, too near 1. Do "
                               "not know what to do.");
      } else if (radius_inter[0] > 1 || radius_inter[0] <= 0) {
         logger(Logger::ERROR, "[generate_complaiant_Grid2D] Radial interval with only one point, outside the "
                               "allowed range of (0,1).");
      } else {
         if (radius_g_size.size() != 1) {
            logger(Logger::ERROR,
                   std::format("[generate_complaiant_Grid2D] Specified many points (or none) for sub-intervals "
                               "(vector size = {:d}), but only one "
                               "interval point is give. Do not know what to do.",
                               radius_g_size.size()));
         } else {
            logger(Logger::WARNING, std::format("[generate_complaiant_Grid2D] Radial interval with only one point, "
                                                "extending it to be [{:+.10e}, 1].",
                                                radius_inter[0]));
            radius_inter.emplace_back(1);
         }
      }
   }

   if (radius_inter.size() - 1 != radius_g_size.size()) {
      logger(Logger::ERROR,
             std::format("[generate_complaiant_Grid2D] Incompatible sizes of interval subdivision and number of points "
                         "for sub-interval: ({:d}, {:d}). The former should be exactly one more than the latter",
                         radius_inter.size(), radius_g_size.size()));
   }
   if (radius_inter[0] <= 0) {
      logger(Logger::ERROR, std::format("[generate_complaiant_Grid2D] Incorrect lower value for radius: {:+.16e}. It "
                                        "should be a value strictly satisfying: 0 < r_min < 1.",
                                        radius_inter[0]));
   }
   if (radius_inter[radius_inter.size() - 1] > 1) {
      logger(Logger::ERROR, std::format("[generate_complaiant_Grid2D] Incorrect upper value for radius: {:+.16e}. It "
                                        "should be a exaclty 1",
                                        radius_inter[radius_inter.size() - 1]));
   }
   if (!is_near(radius_inter[radius_inter.size() - 1], 1.0)) {
      logger(Logger::WARNING,
             std::format(
                 "[generate_complaiant_Grid2D] Incorrect upper value for radius: {:+.16e}. It "
                 "should be a exaclty 1. I will push the 1 and duplicate last entry in the number of points vector.",
                 radius_inter[radius_inter.size() - 1]));
      radius_inter.push_back(1);
      size_t tmp = radius_g_size[radius_g_size.size() - 1];
      radius_g_size.emplace_back(tmp);
   }

   for (size_t i = 1; i < radius_inter.size(); i++) {
      if (radius_inter[i] <= radius_inter[i - 1]) {
         logger(Logger::ERROR, std::format("[generate_complaiant_Grid2D] Incorrectly ordered vector of radial "
                                           "subdivisions. First problematic entry is ({:d}, {:d}), for which r_{:d} = "
                                           "{:+.10e} <= {:.10e} = r_{:d}",
                                           i - 1, i, i, radius_inter[i], radius_inter[i - 1], i - 1));
      }
   }

   // SECTION: Sanity checks on the maps
   // SUBSECTION: Back and forth
   const double r_check = 0.52;
   const double r_bnf   = radius_to_p_space(radius_to_i_space(r_check));
   if (!is_near(r_check, r_bnf)) {
      logger(Logger::ERROR, std::format("The radial map is not correctly coded. From rho(r({:.2f})) I find: {:+.16e}",
                                        r_check, r_bnf));
   }

   const double f_check = 1.27;
   const double f_bnf   = angle_to_p_space(angle_to_i_space(f_check));
   if (!is_near(f_check, f_bnf)) {
      logger(Logger::ERROR,
             std::format("The angular map is not correctly coded. From phi(f(:.2f)) I find: {:+.16e}", f_check, f_bnf));
   }

   // SUBSECTION: To inter space derivatives
   const double dx        = 1.0e-6;
   const double r_der     = radius_to_i_space_der(r_check);
   const double r_num_der = (radius_to_i_space(r_check + dx) - radius_to_i_space(r_check - dx)) / (2.0 * dx);
   if (!is_near(r_der, r_num_der, dx)) {
      logger(Logger::ERROR, std::format("The radial map derivative is not correctly coded. From r=({:.2f}) I find: "
                                        "{:+.16e} instead of {:+.16e}",
                                        r_check, r_der, r_num_der));
   }
   const double f_der     = angle_to_i_space_der(f_check);
   const double f_num_der = (angle_to_i_space(f_check + dx) - angle_to_i_space(f_check - dx)) / (2.0 * dx);
   if (!is_near(f_der, f_num_der, dx)) {
      logger(Logger::ERROR, std::format("The angular map derivative is not correctly coded. From f=({:.2f}) I find: "
                                        "{:+.16e} instead of {:+.16e}",
                                        f_check, f_der, f_num_der));
   }

   // SUBSECTION: To phys space derivatives
   const double r_tmp       = radius_to_i_space_der(r_check);
   const double rho_der     = radius_to_p_space_der(r_tmp);
   const double rho_num_der = (radius_to_p_space(r_tmp + dx) - radius_to_p_space(r_tmp - dx)) / (2.0 * dx);
   if (!is_near(rho_der, rho_num_der, dx)) {
      logger(Logger::ERROR, std::format("The radial map derivative is not correctly coded. From rho'({:.2f}) I find: "
                                        "{:+.16e} instead of {:+.16e}",
                                        r_tmp, rho_der, rho_num_der));
   }
   const double f_tmp       = angle_to_i_space_der(f_check);
   const double phi_der     = angle_to_p_space_der(f_check);
   const double phi_num_der = (angle_to_p_space(f_tmp + dx) - angle_to_p_space(f_tmp - dx)) / (2.0 * dx);
   if (!is_near(phi_der, phi_num_der, dx)) {
      logger(Logger::ERROR, std::format("The angular map derivative is not correctly coded. From phi'({:.2f}) I find: "
                                        "{:+.16e} instead of {:+.16e}",
                                        f_tmp, phi_der, phi_num_der));
   }

   // SECTION: Done with checks, can initialize and return

   SingleDiscretizationInfo info_radius(radius_inter, radius_g_size, false, radius_to_i_space, radius_to_i_space_der,
                                        radius_to_p_space, radius_to_p_space_der);

   const std::vector<size_t> is_a(6, n_pts_for_angle_sector);
   SingleDiscretizationInfo info_angle({0, 1, 2, 3, 4, 5, 6}, is_a, true, angle_to_i_space, angle_to_i_space_der,
                                       angle_to_p_space, angle_to_p_space_der);

   return Grid2D(info_radius, info_angle);
}

//========================================================
//========================================================
//========================================================
//========================================================
//========================================================
//========================================================

Discretization1D::Discretization1D(const Grid &grid, std::function<double(double)> const &function)
    : _grid(grid), size(grid.size), size_li(grid.size_li), _fj(grid.c_size)
{
   // long int index_li = 0;
   // for (size_t k = 0; k < grid.c_size; k++) {
   //    _fj(index_li++) = function(grid._coord[k]);
   // }

   for (size_t k = 0; k < grid.size; k++) {
      _fj(_grid._from_iw_to_ic[k]) = function(grid._coord[_grid._from_iw_to_ic[k]]);
   }
}

double Discretization1D::interpolate_as_weights(double x) const
{

   const double r = _grid._d_info.to_inter_space(x);

   long int index_li = 0;

   double res = 0;

   for (size_t k = 0; k < _grid.size; k++) {
      res += _fj(_grid._from_iw_to_ic[k]) * _grid._weights[k](r);
   }

   return res;
}

} // namespace Honeycomb