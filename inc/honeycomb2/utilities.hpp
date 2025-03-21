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
      Logger() : has_been_initialized(false) {}
      Logger(const Logger &)            = delete;
      Logger &operator=(const Logger &) = delete;

      Logger(Logger &&)            = delete;
      Logger &operator=(Logger &&) = delete;

      static Logger instance;

      enum LEVEL { INFO = 0, WARNING = 1, ERROR = 2 };

      void init(std::string filename, int min_lev = 2, bool to_stderr = true)
      {
         if (to_stderr) fileStream = stderr;
         else fileStream = std::fopen(filename.c_str(), "w");
         min_log_level        = min_lev;
         has_been_initialized = true;
      }

      void operator()(LEVEL l, const std::string &message)
      {
         if (!has_been_initialized) init("EMPTY", 2);
         switch (l) {
         case INFO:
            std::fprintf(fileStream, "\033[0;32m[INFO]\033[0m    ");
            break;
         case WARNING:
            std::fprintf(fileStream, "\033[0;33m[WARNING]\033[0m ");
            break;
         case ERROR:
            std::fprintf(fileStream, "\033[0;31m[ERROR]\033[0m   ");
            break;
         default:
            break;
         }
         std::fprintf(fileStream, "%s\n", message.c_str());
         if (l == ERROR) exit(-1);
      }
};

extern Logger logger;

void to_lower_case(std::string &s);

//_________________________________________________________________________________
std::string to_lower_case(std::string s);

//_________________________________________________________________________________
void setup_logging(std::string outpath, std::string outname, bool to_stderr, size_t min_log_level);
} // namespace Honeycomb
