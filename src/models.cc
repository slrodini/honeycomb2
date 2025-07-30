//
// Author(s): Simone Rodini   <rodini.simone.luigi@gmail.com>
//            Arianna Vercesi <arianna.vercesi01@universitadipavia.it>
//

#include <honeycomb2/solution.hpp>
#include <string>

namespace _private_models
{
using namespace Honeycomb;
InputModel OriginalModel()
{

   // Small utilities
   auto max3 = [](double a, double b, double c) {
      return std::max(a, std::max(b, c));
   };

   auto dome = [](double x1, double x2, double x3) {
      return (1 - x1 * x1) * (1 - x2 * x2) * (1 - x3 * x3);
   };

   Honeycomb::InputModel model;

   model.SetModel(InputModel::T_UP, [=](double x1, double x2, double x3) {
      return dome(x1, x2, x3) * cos(4.0 * x2);
   });

   model.SetModel(InputModel::T_DN, [=](double x1, double x2, double x3) {
      double temp = dome(x1, x2, x3);
      return (2 - cos(3 * M_PI * temp)) * temp;
   });

   model.SetModel(InputModel::T_ST, [=](double x1, double x2, double x3) {
      double temp = dome(x1, x2, x3);
      return -0.3 * (2 - cos(3 * M_PI * temp)) * temp;
   });

   model.SetModel(InputModel::DT_UP, [=](double x1, double x2, double x3) {
      return (sin(x2 * M_PI) + 4 * (x1 * x1 - x3 * x3)) * dome(x1, x2, x3);
   });

   model.SetModel(InputModel::DT_DN, [=](double x1, double x2, double x3) {
      double r = max3(fabs(x1), fabs(x2), fabs(x3));
      return sin(x2 * M_PI) * (2 - 2 * cos(3 * M_PI * dome(x1, x2, x3))) / sqrt(r);
   });

   model.SetModel(InputModel::DT_ST, [=](double x1, double x2, double x3) {
      double r = max3(fabs(x1), fabs(x2), fabs(x3));
      return -0.3 * sin(x2 * M_PI) * (2 - 2 * cos(3 * M_PI * dome(x1, x2, x3))) / sqrt(r);
   });

   model.SetModel(InputModel::T_P_GL, [=](double x1, double x2, double x3) {
      double r = max3(fabs(x1), fabs(x2), fabs(x3));
      return sin(x1 - x3) * dome(x1, x2, x3) * sqrt(r);
   });

   model.SetModel(InputModel::T_M_GL, [=](double x1, double x2, double x3) {
      double r = max3(fabs(x1), fabs(x2), fabs(x3));
      return cos(x1 - x3) * dome(x1, x2, x3) * sqrt(r);
   });

   return model;
}

InputModel LFWA_asymptotic()
{
   InputModel model;

   model.SetModel(InputModel::T_UP, [](double x1, double x2, double x3) -> double {
      double expr1;
      if (x1 >= 0 || x3 <= 0) return 0;

      if (x2 <= 0) {
         expr1 = -1.2533141373155001 * pow(x2, 2) * (-9.89888931329275 * x1 * pow(-1. + x3, 3));
      } else {
         expr1 = -(1.2533141373155001 * pow(x2, 2) * (-9.89888931329275 * pow(1 + x1, 3) * x3));
      }

      return expr1;
   });

   model.SetModel(InputModel::T_DN, [](double x1, double x2, double x3) -> double {
      double expr1;
      if (x1 >= 0 || x3 <= 0) return 0;

      if (x2 <= 0) {
         expr1 = -(1.2533141373155001 * pow(x2, 2) * (9.898889313292752 * x1 * pow(-1. + x3, 3)));
      } else {
         expr1 = -(1.2533141373155001 * pow(x2, 2) * (9.898889313292752 * pow(1 + x1, 3) * x3));
      }

      return expr1;
   });

   model.SetModel(InputModel::DT_UP, [](double x1, double x2, double x3) -> double {
      double expr1;
      if (x1 >= 0 || x3 <= 0) return 0;

      if (x2 <= 0) {
         expr1 = 1.2533141373155001 * pow(x2, 2) * (-12.408466885676829 * x1 * pow(-1. + x3, 3));
      } else {
         expr1 = 1.2533141373155001 * pow(x2, 2) * (+12.408466885676829 * pow(1 + x1, 3) * x3);
      }

      return expr1;
   });

   model.SetModel(InputModel::DT_DN, [](double x1, double x2, double x3) -> double {
      double expr1;
      if (x1 >= 0 || x3 <= 0) return 0;

      if (x2 <= 0) {
         expr1 = 1.2533141373155001 * pow(x2, 2) * (4.879734168524594 * x1 * pow(-1. + x3, 3));
      } else {
         expr1 = 1.2533141373155001 * pow(x2, 2) * (-4.879734168524594 * pow(1 + x1, 3) * x3);
      }

      return expr1;
   });

   return model;
}

InputModel LFWA_fitted()
{
   InputModel model;

   model.SetModel(InputModel::T_UP, [](double x1, double x2, double x3) -> double {
      double expr1;
      if (x1 >= 0 || x3 <= 0) return 0;

      if (x2 <= 0) {
         expr1 = -1.2533141373155001
               * ((-0.2517483540738777
                   * (-0.6347433199714827 * pow(-x1, 0.465) * pow(-x2, 1.403) * pow(1 - x3, 2.213)
                          * pow(x3, 1.065)
                      + 0.8039605205790323 * pow(-x1, 0.299) * pow(-x2, 1.721) * pow(1 - x3, 2.251)
                            * pow(x3, 1.065)))
                      / x3
                  + (0.5034967081477554
                     * (0.40198026028951617 * pow(-x1, 0.299) * pow(-x2, 1.721) * pow(1 - x3, 2.251)
                            * pow(x3, 1.065)
                        - 0.5333399196727042 * pow(-x1, 0.299) * pow(-x2, 1.403)
                              * pow(1 - x3, 2.3790000000000004) * pow(x3, 1.065)
                        + 0.25932516312332793 * pow(-x1, 0.261) * pow(-x2, 1.403)
                              * pow(1 - x3, 2.417) * pow(x3, 1.065)
                        + 1.810493702406733 * pow(-x1, 0.299) * pow(-x2, 1.721)
                              * pow(1 - x3, 1.7279999999999998) * pow(x3, 1.588)
                        - 0.8360701959800317 * pow(-x1, 0.299) * pow(-x2, 1.403)
                              * pow(1 - x3, 1.856) * pow(x3, 1.588)
                        + 0.8051538829358363 * pow(-x1, 0.261) * pow(-x2, 1.403)
                              * pow(1 - x3, 1.894) * pow(x3, 1.588)))
                        / x3);
      } else {
         expr1 = -1.2533141373155001
               * ((-0.5034967081477554
                   * (0.8051538829358363 * pow(-x1, 1.588) * pow(1 + x1, 1.894) * pow(x2, 1.403)
                          * pow(x3, 0.261)
                      + 0.25932516312332793 * pow(-x1, 1.065) * pow(1 + x1, 2.417) * pow(x2, 1.403)
                            * pow(x3, 0.261)
                      - 0.8360701959800317 * pow(-x1, 1.588) * pow(1 + x1, 1.856) * pow(x2, 1.403)
                            * pow(x3, 0.299)
                      - 0.5333399196727042 * pow(-x1, 1.065) * pow(1 + x1, 2.3790000000000004)
                            * pow(x2, 1.403) * pow(x3, 0.299)
                      + 1.810493702406733 * pow(-x1, 1.588) * pow(1 + x1, 1.7279999999999998)
                            * pow(x2, 1.721) * pow(x3, 0.299)
                      + 0.40198026028951617 * pow(-x1, 1.065) * pow(1 + x1, 2.251) * pow(x2, 1.721)
                            * pow(x3, 0.299)))
                      / x1
                  + (0.2517483540738777
                     * (0.8039605205790323 * pow(-x1, 1.065) * pow(1 + x1, 2.251) * pow(x2, 1.721)
                            * pow(x3, 0.299)
                        - 0.6347433199714827 * pow(-x1, 1.065) * pow(1 + x1, 2.213) * pow(x2, 1.403)
                              * pow(x3, 0.465)))
                        / x1);
      }

      return expr1;
   });

   model.SetModel(InputModel::T_DN, [](double x1, double x2, double x3) -> double {
      double expr1;
      if (x1 >= 0 || x3 <= 0) return 0;

      if (x2 <= 0) {
         expr1 = (-1.2533141373155001
                  * (-0.2517483540738777
                         * (-1.280514960367045 * pow(-x1, 0.465) * pow(-x2, 1.403)
                                * pow(1 - x3, 2.213) * pow(x3, 1.065)
                            + 1.6079210411580647 * pow(-x1, 0.299) * pow(-x2, 1.721)
                                  * pow(1 - x3, 2.251) * pow(x3, 1.065))
                     - 0.5034967081477554
                           * (0.40198026028951617 * pow(-x1, 0.299) * pow(-x2, 1.721)
                                  * pow(1 - x3, 2.251) * pow(x3, 1.065)
                              + 0.5333399196727042 * pow(-x1, 0.299) * pow(-x2, 1.403)
                                    * pow(1 - x3, 2.3790000000000004) * pow(x3, 1.065))))
               / x3;
      } else {
         expr1 = -1.2533141373155001
               * ((0.5034967081477554
                   * (0.5333399196727042 * pow(-x1, 1.065) * pow(1 + x1, 2.3790000000000004)
                          * pow(x2, 1.403) * pow(x3, 0.299)
                      + 0.40198026028951617 * pow(-x1, 1.065) * pow(1 + x1, 2.251) * pow(x2, 1.721)
                            * pow(x3, 0.299)))
                      / x1
                  + (0.2517483540738777
                     * (1.6079210411580647 * pow(-x1, 1.065) * pow(1 + x1, 2.251) * pow(x2, 1.721)
                            * pow(x3, 0.299)
                        - 1.280514960367045 * pow(-x1, 1.065) * pow(1 + x1, 2.213) * pow(x2, 1.403)
                              * pow(x3, 0.465)))
                        / x1);
      }

      return expr1;
   });

   model.SetModel(InputModel::DT_UP, [](double x1, double x2, double x3) -> double {
      double expr1;
      if (x1 >= 0 || x3 <= 0) return 0;

      if (x2 <= 0) {
         expr1 = (1.2533141373155001
                  * (0.2517483540738777
                         * (-0.6347433199714827 * pow(-x1, 0.465) * pow(-x2, 1.403)
                                * pow(1 - x3, 2.213) * pow(x3, 1.065)
                            + 0.8039605205790323 * pow(-x1, 0.299) * pow(-x2, 1.721)
                                  * pow(1 - x3, 2.251) * pow(x3, 1.065))
                     + 0.5034967081477554
                           * (0.40198026028951617 * pow(-x1, 0.299) * pow(-x2, 1.721)
                                  * pow(1 - x3, 2.251) * pow(x3, 1.065)
                              - 0.5333399196727042 * pow(-x1, 0.299) * pow(-x2, 1.403)
                                    * pow(1 - x3, 2.3790000000000004) * pow(x3, 1.065)
                              + 0.25932516312332793 * pow(-x1, 0.261) * pow(-x2, 1.403)
                                    * pow(1 - x3, 2.417) * pow(x3, 1.065)
                              + 1.810493702406733 * pow(-x1, 0.299) * pow(-x2, 1.721)
                                    * pow(1 - x3, 1.7279999999999998) * pow(x3, 1.588)
                              - 0.8360701959800317 * pow(-x1, 0.299) * pow(-x2, 1.403)
                                    * pow(1 - x3, 1.856) * pow(x3, 1.588)
                              + 0.8051538829358363 * pow(-x1, 0.261) * pow(-x2, 1.403)
                                    * pow(1 - x3, 1.894) * pow(x3, 1.588))))
               / x3;
      } else {
         expr1 = 1.2533141373155001
               * ((0.5034967081477554
                   * (0.8051538829358363 * pow(-x1, 1.588) * pow(1 + x1, 1.894) * pow(x2, 1.403)
                          * pow(x3, 0.261)
                      + 0.25932516312332793 * pow(-x1, 1.065) * pow(1 + x1, 2.417) * pow(x2, 1.403)
                            * pow(x3, 0.261)
                      - 0.8360701959800317 * pow(-x1, 1.588) * pow(1 + x1, 1.856) * pow(x2, 1.403)
                            * pow(x3, 0.299)
                      - 0.5333399196727042 * pow(-x1, 1.065) * pow(1 + x1, 2.3790000000000004)
                            * pow(x2, 1.403) * pow(x3, 0.299)
                      + 1.810493702406733 * pow(-x1, 1.588) * pow(1 + x1, 1.7279999999999998)
                            * pow(x2, 1.721) * pow(x3, 0.299)
                      + 0.40198026028951617 * pow(-x1, 1.065) * pow(1 + x1, 2.251) * pow(x2, 1.721)
                            * pow(x3, 0.299)))
                      / x1
                  + (0.2517483540738777
                     * (0.8039605205790323 * pow(-x1, 1.065) * pow(1 + x1, 2.251) * pow(x2, 1.721)
                            * pow(x3, 0.299)
                        - 0.6347433199714827 * pow(-x1, 1.065) * pow(1 + x1, 2.213) * pow(x2, 1.403)
                              * pow(x3, 0.465)))
                        / x1);
      }

      return expr1;
   });

   model.SetModel(InputModel::DT_DN, [](double x1, double x2, double x3) -> double {
      double expr1;
      if (x1 >= 0 || x3 <= 0) return 0;

      if (x2 <= 0) {
         expr1 = 1.2533141373155001
               * ((0.2517483540738777
                   * (-1.280514960367045 * pow(-x1, 0.465) * pow(-x2, 1.403) * pow(1 - x3, 2.213)
                          * pow(x3, 1.065)
                      + 1.6079210411580647 * pow(-x1, 0.299) * pow(-x2, 1.721) * pow(1 - x3, 2.251)
                            * pow(x3, 1.065)))
                      / x3
                  - (0.5034967081477554
                     * (0.40198026028951617 * pow(-x1, 0.299) * pow(-x2, 1.721) * pow(1 - x3, 2.251)
                            * pow(x3, 1.065)
                        + 0.5333399196727042 * pow(-x1, 0.299) * pow(-x2, 1.403)
                              * pow(1 - x3, 2.3790000000000004) * pow(x3, 1.065)))
                        / x3);
      } else {
         expr1 = 1.2533141373155001
               * ((-0.5034967081477554
                   * (0.5333399196727042 * pow(-x1, 1.065) * pow(1 + x1, 2.3790000000000004)
                          * pow(x2, 1.403) * pow(x3, 0.299)
                      + 0.40198026028951617 * pow(-x1, 1.065) * pow(1 + x1, 2.251) * pow(x2, 1.721)
                            * pow(x3, 0.299)))
                      / x1
                  + (0.2517483540738777
                     * (1.6079210411580647 * pow(-x1, 1.065) * pow(1 + x1, 2.251) * pow(x2, 1.721)
                            * pow(x3, 0.299)
                        - 1.280514960367045 * pow(-x1, 1.065) * pow(1 + x1, 2.213) * pow(x2, 1.403)
                              * pow(x3, 0.465)))
                        / x1);
      }

      return expr1;
   });

   return model;
}

// =============================================================================
InputModel SiversExtension()
{

   // Small utilities
   auto max3 = [](double a, double b, double c) {
      return std::max(a, std::max(b, c));
   };

   auto dome = [](double x1, double x2, double x3) {
      return (1 - x1 * x1) * (1 - x2 * x2) * (1 - x3 * x3) / pow(1.0 - pow(x1 - x3, 2) / 4., 2.);
   };

   auto angle = [=](double x1, double x2, double x3) -> double {
      double r = max3(std::fabs(x1), std::fabs(x2), std::fabs(x3));
      if (x1 > 0 && x2 >= 0 && x3 < 0) return x2 / r;
      if (x1 <= 0 && x2 > 0 && x3 < 0) return 1. - x1 / r;
      if (x1 < 0 && x2 > 0 && x3 >= 0) return 3. - x2 / r;
      if (x1 < 0 && x2 <= 0 && x3 > 0) return 3. - x2 / r;
      if (x1 >= 0 && x2 < 0 && x3 > 0) return 4. + x1 / r;
      if (x1 > 0 && x2 < 0 && x3 <= 0) return 6. + x2 / r;
      return 0.;
   };

   struct LocParam {
      double Nsiv_a, G_a, alpha_a, beta_a, Atilde_a, Btilde_a;
   };

   LocParam tmd_u = {
       .Nsiv_a   = 0.42,
       .G_a      = 0.0199,
       .alpha_a  = 0.41,
       .beta_a   = 1.66,
       .Atilde_a = -0.58,
       .Btilde_a = 1.12,
   };
   LocParam tmd_d = {
       .Nsiv_a   = -1.0,
       .G_a      = 0.000558,
       .alpha_a  = 0.94,
       .beta_a   = 10.00,
       .Atilde_a = -0.78,
       .Btilde_a = 0.98,
   };

   LocParam tmd_s = {
       .Nsiv_a   = 0.28,
       .G_a      = 0.00242,
       .alpha_a  = 0.61,
       .beta_a   = 8.23,
       .Atilde_a = -1.44,
       .Btilde_a = 0.92,
   };

   struct PDFParam {
      double N_a, a_a, b_a, A_a, B_a;
   };
   PDFParam u_valence = {
       .N_a = 0.5889,
       .a_a = 0.3444,
       .b_a = 3.7312,
       .A_a = -0.1740,
       .B_a = 17.997,
   };
   PDFParam d_valence = {
       .N_a = 0.2585,
       .a_a = 0.2951,
       .b_a = 4.8682,
       .A_a = -1.0552,
       .B_a = 26.536,
   };

   PDFParam db_m_ub = {
       .N_a = 7.2847,
       .a_a = 1.2773,
       .b_a = 18.756,
       .A_a = -6.3187,
       .B_a = 18.306,
   };

   PDFParam db_p_ub = {
       .N_a = 0.2295,
       .a_a = -0.1573,
       .b_a = 8.8819,
       .A_a = 0.8704,
       .B_a = 8.2179,
   };

   auto pdf = [](double x, PDFParam p) {
      return p.N_a * pow(x, p.a_a - 1) * pow(1.0 - x, p.b_a) * (1.0 + p.A_a * sqrt(x) + p.B_a);
   };

   auto pdf_u = [=](double x) {
      return pdf(x, u_valence) + 0.5 * (pdf(x, db_p_ub) - pdf(x, db_m_ub));
   };
   auto pdf_ub = [=](double x) {
      return -pdf(x, u_valence) + pdf_u(x);
   };

   auto pdf_d = [=](double x) {
      return pdf(x, d_valence) + 0.5 * (pdf(x, db_p_ub) + pdf(x, db_m_ub));
   };
   auto pdf_db = [=](double x) {
      return -pdf(x, d_valence) + pdf_d(x);
   };

   auto pdf_s = [=](double x) {
      return 0.5 * pdf(x, db_p_ub);
   };
   auto pdf_sb = [=](double x) {
      return 0.5 * pdf(x, db_p_ub);
   };

   auto Kx = [](double x) {
      double N1     = 0.285;
      double k1     = 2.98;
      double k2     = 0.173;
      double lambda = 0.38;
      double Msq    = pow(0.938, 2);
      double xhat   = 0.1;
      double g1     = N1 * pow(x / xhat, k2) * pow((1 - x) / (1 - xhat), k1);

      double M1sq    = 0.44 * 0.44;
      double lambdaS = 2.;

      double den
          = 2. * M_PI * M_PI * Msq * (1. + lambda * g1) * pow(g1 + M1sq, 2) * (1. + lambdaS * M1sq);
      double num = g1 * M1sq
                 * (1. + 2. * (lambda + lambdaS) * g1 * M1sq / (g1 + M1sq)
                    + 6 * lambda * lambdaS * pow(g1 * M1sq / (g1 + M1sq), 2.));

      return num / den;
   };

   auto Ta = [=](double x, LocParam p, std::function<double(double)> _pdf) {
      double t1 = p.Nsiv_a * Kx(x) * pow(x, p.alpha_a) * pow(1.0 - x, p.beta_a) / p.G_a;
      double t2 = (1. - p.Btilde_a + p.Atilde_a * x + 2. * p.Btilde_a * x * x);
      double fa = _pdf(x);
      return t1 * t2 * fa;
   };

   auto f_up = [=](double x) {
      return Ta(x, tmd_u, pdf_u);
   };
   auto g_up = [=](double x) {
      return Ta(x, tmd_u, pdf_ub);
   };

   auto f_dn = [=](double x) {
      return Ta(x, tmd_d, pdf_d);
   };
   auto g_dn = [=](double x) {
      return Ta(x, tmd_d, pdf_db);
   };

   auto f_st = [=](double x) {
      return Ta(x, tmd_s, pdf_s);
   };
   auto g_st = [=](double x) {
      return Ta(x, tmd_s, pdf_sb);
   };

   auto extension = [=](double x1, double x2, double x3, std::function<double(double)> _f,
                        std::function<double(double)> _g) -> double {
      double pref = dome(x1, x2, x3);
      double phi  = angle(x1, x2, x3);
      if (0 <= phi && phi < 1) return pref * (1. - phi) * _g(0.5 * (x1 - x3));
      if (1 <= phi && phi < 2) return 0.;
      if (2 <= phi && phi < 3) return pref * (phi - 2.) * _f(0.5 * (x3 - x1));
      if (3 <= phi && phi < 4) return pref * (4. - phi) * _f(0.5 * (x3 - x1));
      if (4 <= phi && phi < 5) return 0.;
      if (5 <= phi && phi <= 6) return pref * (phi - 5.) * _g(0.5 * (x1 - x3));
      return 0.;
   };

   Honeycomb::InputModel model;

   model.SetModel(InputModel::T_UP, [=](double x1, double x2, double x3) {
      double r = max3(std::fabs(x1), std::fabs(x2), std::fabs(x3));
      if (std::fabs(1. - r) < 1.0e-15) return 0.;
      return extension(x1, x2, x3, f_up, g_up);
   });

   model.SetModel(InputModel::T_DN, [=](double x1, double x2, double x3) {
      double r = max3(std::fabs(x1), std::fabs(x2), std::fabs(x3));
      if (std::fabs(1. - r) < 1.0e-15) return 0.;
      return extension(x1, x2, x3, f_dn, g_dn);
   });

   model.SetModel(InputModel::T_ST, [=](double x1, double x2, double x3) {
      double r = max3(std::fabs(x1), std::fabs(x2), std::fabs(x3));
      if (std::fabs(1. - r) < 1.0e-15) return 0.;
      return extension(x1, x2, x3, f_st, g_st);
   });

   return model;
}

} // namespace _private_models

namespace Honeycomb
{

std::map<std::string, std::function<InputModel()>> PreImplementedModels::_models
    = {{"pim_original", _private_models::OriginalModel},
       {"pim_LFWA_asymptotic", _private_models::LFWA_asymptotic},
       {"pim_LFWA_fitted", _private_models::LFWA_fitted},
       {"pim_SiversExtension", _private_models::SiversExtension}};

InputModel PreImplementedModels::GetModel(std::string key)
{
   if (_models.find(key) != _models.end()) {
      return _models.at(key)();
   } else {
      logger(Logger::ERROR, "Model: <" + key + "> is not available.");
      return InputModel();
   }
}

bool PreImplementedModels::QueryModelAvailability(std::string model)
{
   return _models.find(model) != _models.end();
}

void PreImplementedModels::AddModel(std::string name, std::function<InputModel()> model)
{
   if (_models.find(name) == _models.end()) {
      _models[name] = model;
   } else {
      size_t i = 2;
      std::string tmp_name;
      while (true) {
         tmp_name = name + "_" + std::to_string(i);
         if (_models.find(tmp_name) == _models.end()) {
            _models[tmp_name] = model;
         } else {
            i++;
         }
      }
      logger(Logger::WARNING, "Model: <" + name + "> is already present. "
                                  + "I will add this model with key: " + tmp_name);
   }
}

} // namespace Honeycomb