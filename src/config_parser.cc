#include <honeycomb2/config_parser.hpp>
#include <honeycomb2/utilities.hpp>

namespace Honeycomb
{

std::string read_file_to_str(const std::string &file_path)
{
   std::FILE *fp = std::fopen(file_path.c_str(), "rb");
   if (!fp)
      logger(Logger::ERROR,
             std::format("ConfigParser::_read_entire_file Impossible to open the file <{:s}>",
                         file_path));

   std::fseek(fp, 0L, SEEK_END);
   long file_size = std::ftell(fp);
   std::rewind(fp);

   std::string content(file_size + 1, '\0');

   size_t n_char = std::fread(content.data(), 1, file_size, fp);
   if (n_char != (size_t)file_size)
      logger(Logger::ERROR, std::format("ConfigParser::_read_entire_file Read {:d} character where "
                                        "size of file <{:s}> is {:d}",
                                        n_char, file_path, file_size));
   std::fclose(fp);

   return content;
}

std::vector<std::string> ConfigParser::_get_lines(const std::string &str)
{
   if (str.size() == 0) return {};

   bool skipping = true;

   std::vector<std::string> tokens;
   tokens.emplace_back("");

   for (size_t i = 0; i < str.size(); i++) {
      bool split = (str[i] == '\n');

      if (split && skipping) continue;
      if (split && !skipping) {
         tokens.emplace_back("");
         skipping = true;
      }
      if (!split) {
         tokens.back() += str[i];
         skipping       = false;
      }
   }

   return tokens;
}

std::vector<std::string> ConfigParser::_tokenize_string(const std::string &str,
                                                        const std::string &delim)
{
   if (str.size() == 0) return {};

   bool skipping = true;

   std::vector<std::string> tokens;
   tokens.emplace_back("");

   for (size_t i = 0; i < str.size(); i++) {
      if (str[i] == '\"') {
         i++;
         while (i < str.size() && str[i] != '\"') {
            tokens.back() += str[i];
            i++;
         }
         if (i == str.size()) {
            logger(Logger::ERROR, "ConfigParser::_tokenize_string unmatched \" in config file.");
         } else if (i == str.size() - 1) {
            break;
         } else {
            i++;
         }
      }
      if (i >= str.size()) break;

      bool split = false;
      for (size_t j = 0; j < delim.size(); j++) {
         split = split || (str[i] == delim[j]);
         if (split) break;
      }

      if (split && skipping) continue;
      if (split && !skipping) {
         tokens.emplace_back("");
         skipping = true;
      }
      if (!split) {
         tokens.back() += str[i];
         skipping       = false;
      }
   }

   return tokens;
}

ConfigParser::ConfigParser(const std::string &content)
{
   std::vector<std::string> lines = _get_lines(content);

   std::string delim(" ,;:\n\t=");

   for (auto &line : lines) {
      if (_check_skip_line(line) != 0) continue;
      _clean_trailing_whitespace(line);

      for (size_t j = 0; j < line.size(); j++) {
         if (line[j] == '#') {
            line.resize(j);
            break;
         }
      }

      {
         std::vector<std::string> tokens = _tokenize_string(line, delim);
         std::string key                 = tokens[0];
         std::vector<std::string> values(tokens.begin() + 1, tokens.end());
         if (false == _key_does_not_exist(key)) {
            logger(Logger::WARNING, "ConfigParser::ConfigParser Key " + key
                                        + "already exists. I will overwrite the values.");
         }
         _options[key] = values;
      }
   }
}

static inline bool check_not_num(char c)
{
   return c < '0' || c > '9';
}

int8_t ConfigParser::_check_number_nonsense(const std::string &tk, InputType_e fmt)
{
   switch (fmt) {
   case INTEGER: {
      if (check_not_num(tk[0]) && tk[0] != '+' && tk[0] != '-') {
         return 1;
      }
      for (size_t i = 1; i < tk.size(); i++) {
         if (tk[i] == '\0') continue;
         if (check_not_num(tk[i])) {
            return 1;
         }
      }
      break;
   }
   case U_INTEGER: {
      if (check_not_num(tk[0]) && tk[0] != '+') {
         return 1;
      }
      for (size_t i = 1; i < tk.size(); i++) {
         if (tk[i] == '\0') continue;
         if (check_not_num(tk[i])) {
            return 1;
         }
      }
      break;
   }
   case FLOAT: {

      bool dot = 0;
      if (tk[0] && tk[0] != '+' && tk[0] != '-') {
         if (tk[0] != '.') {
            return 1;
         } else {
            dot = true;
         }
      }
      for (size_t i = 1; i < tk.size(); i++) {
         if (tk[i] == '\0') continue;
         if (check_not_num(tk[i])) {
            if (tk[i] != '.') {
               return 1;
            } else if (!dot) {
               dot = true;
            } else {
               return 1;
            }
         }
      }
      break;
   }
   case EXPONENTIAL: {
      size_t l = tk.size();
      bool dot = false;
      if (check_not_num(tk[0]) && tk[0] != '+' && tk[0] != '-') {
         if (tk[0] != '.') {
            return 1;
         } else {
            dot = true;
         }
      }
      bool e    = false;
      bool sign = false;
      for (size_t i = 1; i < l; i++) {
         if (tk[i] == '\0') continue;
         if (check_not_num(tk[i])) {
            if (tk[i] != '.' && tk[i] != 'e' && tk[i] != 'E' && tk[i] != '+' && tk[i] != '-') {
               return 1;
            } else {
               switch (tk[i]) {
               case '.': {
                  if (!dot) {
                     dot = true;
                     break;
                  } else goto error_goto;
                  break;
               }
               case '+': {
                  if (!sign) {
                     sign = true;
                     break;
                  } else goto error_goto;
                  break;
               }
               case '-': {
                  if (!sign) {
                     sign = true;
                     break;
                  } else goto error_goto;
                  break;
               }
               case 'e':
               case 'E': {
                  if (!e) {
                     e = true;
                     break;
                  } else goto error_goto;
                  break;
               }

               default: {
               error_goto:
                  return 1;
               }
               }
            }
         }
      }
      break;
   }
   default:
      break;
   }
   return 0;
}

bool ConfigParser::_key_does_not_exist(const std::string &key)
{
   return _options.find(key) == _options.end();
}

void ConfigParser::_clean_trailing_whitespace(std::string &line)
{
   for (long int i = static_cast<long int>(line.size()) - 1; i >= 0; i--) {
      if (line[i] == ' ' || line[i] == '\t' || line[i] == '\n' || line[i] == '\0') {
         line.erase(i, 1);
      } else {
         return;
      }
   }
}

int8_t ConfigParser::_check_skip_line(std::string &str)
{
   for (long int i = 0; i < static_cast<long int>(str.size()); i++) {
      if (str[i] == '\n' || str[i] == '#' || str[i] == '\0') return 2;
      if (str[i] == ' ' || str[i] == '\t') {
         str.erase(i, 1);
         i--;
      } else {
         return 0;
      }
   }
   return 1;
}

} // namespace Honeycomb