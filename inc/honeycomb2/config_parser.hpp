#ifndef HC2_CONFIG_PARSER_HPP
#define HC2_CONFIG_PARSER_HPP

#include <honeycomb2/default.hpp>
#include <honeycomb2/utilities.hpp>

namespace Honeycomb
{

std::string read_file_to_str(const std::string &file_path);

namespace custom_concepts
{

template <typename T>
concept arithmetic_vector = requires {
   typename T::value_type;
} && (std::integral<typename T::value_type> || std::floating_point<typename T::value_type>);

} // namespace custom_concepts

class ConfigParser
{
public:
   ConfigParser(const std::string &content);

   const std::map<std::string, std::vector<std::string>> &GetMap()
   {
      return _options;
   }

   template <std::signed_integral T>
   T GetValue(const std::string &key, size_t entry = 0)
   {
      if (_key_does_not_exist(key)) {
         std::cerr << "[ERROR]: Key " << key << "does not exist.";
         return 0;
      }
      if (_check_number_nonsense(_options[key][entry], INTEGER) == 0) {
         T res    = 0;
         auto tmp = std::from_chars(_options[key][entry].data(),
                                    _options[key][entry].data() + _options[key][entry].size(), res);
         if (tmp.ec != std::errc()) {
            std::cerr << "[ERROR]: Parsing failed.\n";
            res = 0;
         }
         return res;
      } else {
         std::cerr << "[ERROR]: Parsing failed, wrong format.\n";
         return 0;
      }
   }

   template <std::unsigned_integral T>
   T GetValue(const std::string &key, size_t entry = 0)
   {
      if (_key_does_not_exist(key)) {
         std::cerr << "[ERROR]: Key " << key << "does not exist.";
         return 0;
      }
      if (_check_number_nonsense(_options[key][entry], U_INTEGER) == 0) {
         T res    = 0;
         auto tmp = std::from_chars(_options[key][entry].data(),
                                    _options[key][entry].data() + _options[key][entry].size(), res);
         if (tmp.ec != std::errc()) {
            std::cerr << "[ERROR]: Parsing failed.\n";
            res = 0;
         }
         return res;
      } else {
         std::cerr << "[ERROR]: Parsing failed, wrong format.\n";
         return 0;
      }
   }

   template <std::floating_point T>
   T GetValue(const std::string &key, size_t entry = 0)
   {
      if (_key_does_not_exist(key)) {
         std::cerr << "[ERROR]: Key " << key << "does not exist.";
         return static_cast<T>(NAN);
      }
      if (_check_number_nonsense(_options[key][entry], EXPONENTIAL) == 0
          || _check_number_nonsense(_options[key][entry], FLOAT) == 0) {
         T res    = 0;
         auto tmp = std::from_chars(_options[key][entry].data(),
                                    _options[key][entry].data() + _options[key][entry].size(), res);
         if (tmp.ec != std::errc()) {
            std::cerr << "[ERROR]: Parsing failed.\n";
            res = static_cast<T>(NAN);
         }
         return res;
      } else {
         std::cerr << "[ERROR]: Parsing failed, wrong format.\n";
         return static_cast<T>(NAN);
      }
   }

   template <custom_concepts::arithmetic_vector V>
   V GetValue(const std::string &key)
   {

      using T = typename V::value_type;

      if (_key_does_not_exist(key)) {
         return {};
      }
      size_t n = _options[key].size();
      V res;
      for (size_t i = 0; i < n; i++) {
         res.emplace_back(GetValue<T>(key, i));
      }
      return res;
   }

   std::string GetValue(const std::string &key, size_t entry = 0)
   {
      if (_key_does_not_exist(key)) {
         logger(Logger::ERROR, "ConfigParser::GetValue Key " + key + " does not exists.");
         return "";
      }
      return _options[key][entry];
   }

   // template <typename T>
   // T GetValue(const std::string &key, size_t entry = 0)
   // {
   //    std::cerr << "[ERROR]: Unknown conversion.\n";
   //    exit(-1);
   // }

private:
   enum InputType_e { INTEGER, U_INTEGER, FLOAT, EXPONENTIAL, STRING };

   std::map<std::string, std::vector<std::string>> _options;

   std::vector<std::string> _tokenize_string(const std::string &str, const std::string &delim);
   std::vector<std::string> _get_lines(const std::string &str);

   void _clean_trailing_whitespace(std::string &line);
   int8_t _check_skip_line(std::string &str);
   int8_t _check_number_nonsense(const std::string &tk, InputType_e fmt);

   bool _key_does_not_exist(const std::string &key);
};

} // namespace Honeycomb
#endif // HC2_CONFIG_PARSER_HPP