#include <honeycomb2/default.hpp>
#include <honeycomb2/utilities.hpp>

namespace Honeycomb
{

Logger logger;

void to_lower_case(std::string &s)
{
   std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) {
      return std::tolower(c);
   });
}

//_________________________________________________________________________________
std::string to_lower_case(std::string s)
{
   std::string ret;
   for (size_t i = 0; i < s.length(); i++)
      ret += std::tolower(s[i]);
   return ret;
}

//_________________________________________________________________________________
void setup_logging(std::string outpath, std::string outname, bool to_stderr, size_t min_log_level)
{

   if (!std::filesystem::exists(outpath)) std::filesystem::create_directories(outpath);
   else {
      if (std::filesystem::is_regular_file(outpath)) {
         std::filesystem::rename(outpath, outpath + ".bak");
         std::filesystem::create_directories(outpath);
      }
   }
   logger.init(outpath + "/" + outname + "_logger.txt", 0);
}
} // namespace Honeycomb
