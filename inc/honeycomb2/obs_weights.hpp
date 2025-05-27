#ifndef G2_WEIGHTS_HPP
#define G2_WEIGHTS_HPP

#include <honeycomb2/discretization.hpp>
#include <honeycomb2/solution.hpp>
#include <honeycomb2/gauss_kronrod.hpp>

namespace Honeycomb
{
struct G2Weights {

   G2Weights(double _xBj, const Grid2D &_grid, bool to_compute = true, double int_e_r = 1.0e-8,
             double int_e_a = 1.0e-8);

   static Eigen::VectorXd GetWeights(double _xBj, const Grid2D &_grid, double int_e_r = 1.0e-8,
                                     double int_e_a = 1.0e-8);

   double xBj;
   const Grid2D &grid;
   Eigen::VectorXd weights;

   template <class Archive>
   void save(Archive &archive) const
   {
      archive(grid);
      archive(xBj);
      archive(weights.size());
      for (long int i = 0; i < weights.size(); i++)
         archive(weights[i]);
   }

   template <class Archive>
   void load(Archive &archive)
   {
      Grid2D tmp_grid;
      archive(tmp_grid);
      check_grid_compatibility(tmp_grid, grid);

      archive(xBj);
      long int s;
      archive(s);
      weights = Eigen::VectorXd::Zero(s);
      for (long int i = 0; i < s; i++)
         archive(weights[i]);
   }
};

// Weights for the d_2 observable, defined as
// 3\int_0^1 dx x^2 g_2(x, Q) = \int Dx \mathfrak{S}^+
struct D2Weights {

   D2Weights(const Grid2D &_grid, bool to_compute = true, double int_e_r = 1.0e-8, double int_e_a = 1.0e-8);

   void GetWeights(double int_e_r = 1.0e-8, double int_e_a = 1.0e-8);

   double ComputeSingleQuark(const Eigen::VectorXd &_f) const;

   const Grid2D &grid;
   Eigen::VectorXd weights;

   template <class Archive>
   void save(Archive &archive) const
   {
      archive(grid);
      archive(weights.size());
      for (long int i = 0; i < weights.size(); i++)
         archive(weights[i]);
   }

   template <class Archive>
   void load(Archive &archive)
   {
      Grid2D tmp_grid;
      archive(tmp_grid);
      check_grid_compatibility(tmp_grid, grid);

      long int s;
      archive(s);
      weights = Eigen::VectorXd::Zero(s);
      for (long int i = 0; i < s; i++)
         archive(weights[i]);
   }
};

struct D2WeightsCutted {

   D2WeightsCutted(const Grid2D &_grid, bool to_compute = true, double int_e_r = 1.0e-8,
                   double int_e_a = 1.0e-8);

   void GetWeights(double int_e_r = 1.0e-8, double int_e_a = 1.0e-8);

   double ComputeSingleQuark(const Eigen::VectorXd &_f) const;
   double ComputeSingleQuark_NoCorrections(const Eigen::VectorXd &_f) const;

   const Grid2D &grid;
   Eigen::VectorXd weights;
   double center_approx;

   template <class Archive>
   void save(Archive &archive) const
   {
      archive(grid);
      archive(center_approx);
      archive(weights.size());
      for (long int i = 0; i < weights.size(); i++)
         archive(weights[i]);
   }

   template <class Archive>
   void load(Archive &archive)
   {
      Grid2D tmp_grid;
      archive(tmp_grid);
      check_grid_compatibility(tmp_grid, grid);
      archive(center_approx);

      long int s;
      archive(s);
      weights = Eigen::VectorXd::Zero(s);
      for (long int i = 0; i < s; i++)
         archive(weights[i]);
   }
};

// Weights for computing d2 as a cutted integral: \int_a^b dx x^2 g2(x, Q)
struct D2WeightsPartialIntegral {

   D2WeightsPartialIntegral(const Grid2D &_grid, double a, double b, bool to_compute = true,
                            double int_e_r = 1.0e-8, double int_e_a = 1.0e-8);

   void GetWeights(double int_e_r = 1.0e-8, double int_e_a = 1.0e-8);

   double ComputeSingleQuark(const Eigen::VectorXd &_f) const;

   const Grid2D &grid;
   const double a;
   const double b;
   Eigen::VectorXd weights;

   template <class Archive>
   void save(Archive &archive) const
   {
      archive(grid);
      archive(a, b);
      archive(weights.size());
      for (long int i = 0; i < weights.size(); i++)
         archive(weights[i]);
   }

   template <class Archive>
   void load(Archive &archive)
   {
      Grid2D tmp_grid;
      archive(tmp_grid);
      check_grid_compatibility(tmp_grid, grid);
      archive(a, b);

      long int s;
      archive(s);
      weights = Eigen::VectorXd::Zero(s);
      for (long int i = 0; i < s; i++)
         archive(weights[i]);
   }
};

// Weights for  Efremov-Leader-Teryaev sum rule, defined as
// 2\int_0^1 dx x g_2(x, Q) = \int_0^1 dx_2 \int_{-x_2}^0 dx_1 (-x_1/x_2^2)\mathfrak{S}^+
struct ELTWeights {

   ELTWeights(const Grid2D &_grid, bool to_compute = true, double int_e_r = 1.0e-8, double int_e_a = 1.0e-8);

   void GetWeights(double int_e_r = 1.0e-8, double int_e_a = 1.0e-8);

   double ComputeSingleQuark(const Eigen::VectorXd &_f) const;
   double ComputeSingleQuark_NoCorrections(const Eigen::VectorXd &_f) const;

   const Grid2D &grid;
   Eigen::VectorXd weights;
   double center_approx;

   template <class Archive>
   void save(Archive &archive) const
   {
      archive(grid);
      archive(center_approx);
      archive(weights.size());
      for (long int i = 0; i < weights.size(); i++)
         archive(weights[i]);
   }

   template <class Archive>
   void load(Archive &archive)
   {
      Grid2D tmp_grid;
      archive(tmp_grid);
      check_grid_compatibility(tmp_grid, grid);
      archive(center_approx);

      long int s;
      archive(s);
      weights = Eigen::VectorXd::Zero(s);
      for (long int i = 0; i < s; i++)
         archive(weights[i]);
   }
};

} // namespace Honeycomb

#endif
