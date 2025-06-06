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
} // namespace _private_models

namespace Honeycomb
{

std::map<std::string, std::function<InputModel()>> PreImplementedModels::_models
    = {{"pim_original", _private_models::OriginalModel}};

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