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

} // namespace _private_models

namespace Honeycomb
{

std::map<std::string, std::function<InputModel()>> PreImplementedModels::_models
    = {{"pim_original", _private_models::OriginalModel},
       {"pim_LFWA_asymptotic", _private_models::LFWA_asymptotic},
       {"pim_LFWA_fitted", _private_models::LFWA_fitted}};

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