#include <honeycomb2/honeycomb2.hpp>

int main()
{

   double fnc_elapsed           = 0;
   Honeycomb::timer::mark begin = Honeycomb::timer::now();
   Honeycomb::timer::mark end   = Honeycomb::timer::now();

   // fnc_elapsed                = Honeycomb::timer::elapsed_ns(end, begin);

   // const double Nc = 3; // NC = 1 for tests

   const size_t n    = 8;
   const double rmin = 0.001;
   Honeycomb::Grid2D grid
       = Honeycomb::generate_compliant_Grid2D(n, {rmin, 0.1, 0.4, 1}, {12, 8, 7});

   Honeycomb::logger(Honeycomb::Logger::INFO,
                     std::format("Total grid size: {:d}x{:d}", grid.size, grid.size));
   Honeycomb::logger(Honeycomb::Logger::INFO,
                     std::format("Total x123 size: {:d}x{:d}", grid.c_size, grid.c_size));

   Honeycomb::Discretization discr(grid);
   Honeycomb::InputModel model = Honeycomb::PreImplementedModels::GetModel("pim_original");
   Honeycomb::Solution sol(&discr, model, 3);

   sol.RotateToPhysicalBasis();

   begin = Honeycomb::timer::now();
   Honeycomb::VectorG2Weights g2_weights(grid, true, 1.0e-10);

   end         = Honeycomb::timer::now();
   fnc_elapsed = Honeycomb::timer::elapsed_ms(end, begin);
   std::cout << "Elapsed (ms): " << fnc_elapsed << std::endl;

   std::vector<double> xBjs = {0.005, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9};
   std::vector<double> exact
       = {5.817094557526623,   1.4640797528367222,  0.22241967743949864, -0.7828228062266303,
          -1.0343438687806414, -0.8910383648576137, -0.5227800633836773, -0.07281606455573503,
          0.3165268263751741,  0.5099529124925853,  0.4052714626463452};

   Eigen::VectorXd F_test_1 = sol._distr_p[2];

   for (size_t j = 0; j < xBjs.size(); j++) {
      double res_g2 = exact[j];
      double tmp    = g2_weights.interpolate(xBjs[j]).dot(F_test_1);

      if (!Honeycomb::is_near(res_g2, tmp, 1.03 - 3)) return 1;
   }

   return 0;
}