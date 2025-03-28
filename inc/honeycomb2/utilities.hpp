#pragma once

#include <honeycomb2/default.hpp>

namespace Honeycomb
{

struct Logger {
private:
   // Private constructor prevents external instantiation
   std::FILE *fileStream;
   bool has_been_initialized;
   int min_log_level = 2;

public:
   Logger() : has_been_initialized(false)
   {
   }
   Logger(const Logger &)            = delete;
   Logger &operator=(const Logger &) = delete;

   Logger(Logger &&)            = delete;
   Logger &operator=(Logger &&) = delete;

   static Logger instance;

   enum LEVEL { INFO = 0, WARNING = 1, ERROR = 2 };
   enum PREV_STATE { UNDEF = -1, TO_FILE = 0, TO_STD = 1 };

   // Calling init multiple times with the same filepath closes the old file AND OVERWRITES IT.
   // Ideally `init` is called exactly once, via the setup_logging utility at the beginning of the program
   void init(std::string filename, int min_lev = 2, bool to_stderr = true)
   {
      if (state == PREV_STATE::UNDEF) {
         if (to_stderr) fileStream = stderr;
         else fileStream = std::fopen(filename.c_str(), "w");
         has_been_initialized = true;
      } else {
         if (state == TO_FILE) std::fclose(fileStream);
         if (to_stderr) fileStream = stderr;
         else fileStream = std::fopen(filename.c_str(), "w");
      }
      min_log_level = min_lev;
      state         = to_stderr ? PREV_STATE::TO_STD : PREV_STATE::TO_FILE;
   }

   void operator()(LEVEL l, const std::string &message)
   {
      if (!has_been_initialized) init("EMPTY", 2);
      switch (l) {
      case LEVEL::INFO:
         if (min_log_level >= LEVEL::INFO) std::fprintf(fileStream, "\033[0;32m[INFO]\033[0m    ");
         break;
      case LEVEL::WARNING:
         if (min_log_level >= LEVEL::WARNING) std::fprintf(fileStream, "\033[0;33m[WARNING]\033[0m ");
         break;
      case LEVEL::ERROR:
         if (min_log_level >= LEVEL::ERROR) std::fprintf(fileStream, "\033[0;31m[ERROR]\033[0m   ");
         break;
      default:
         break;
      }
      std::fprintf(fileStream, "%s\n", message.c_str());
      if (l == LEVEL::ERROR) exit(-1);
   }

private:
   PREV_STATE state = PREV_STATE::UNDEF;
};

extern Logger logger;

void to_lower_case(std::string &s);

//_________________________________________________________________________________
std::string to_lower_case(std::string s);

//_________________________________________________________________________________
void setup_logging(std::string outpath, std::string outname, bool to_stderr, size_t min_log_level);

//_________________________________________________________________________________
inline bool is_near(double a, double b, double eps = 1.0e-15)
{
   return std::fabs(a - b) < eps;
}
} // namespace Honeycomb
