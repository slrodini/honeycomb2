#include <honeycomb2/checksum.hpp>

namespace Honeycomb
{

// CRC-64/WE, without final xor
// Resources: https://www.wolfgang-ehrhardt.de/ ;
// https://reveng.sourceforge.io/crc-catalogue/all.htm
#define POLYNOMIAL 0x42F0E1EBA9EA3693UL
#define MSB_CRC64 0x8000000000000000UL
#define INITIAL_CRC64 0xFFFFFFFFFFFFFFFFUL

//_________________________________________________________________________________
uint64_t crc64(const uint8_t *data, size_t length)
{
   uint64_t crc = INITIAL_CRC64;

   for (size_t j = 0; j < length; j++) {
      uint8_t byte  = data[j];
      crc          ^= static_cast<uint64_t>(byte) << 56;
      for (int i = 0; i < 8; ++i)
         crc = (crc & (MSB_CRC64)) ? (crc << 1) ^ (POLYNOMIAL) : crc << 1;
   }

   // return crc ^ INITIAL_CRC64;
   return crc;
}

//_________________________________________________________________________________
uint64_t compute_file_checksum(const char *file_path)
{
   std::FILE *fp = std::fopen(file_path, "rb");
   if (!fp)
      logger(Logger::ERROR,
             "[compute_file_checksum]: Impossible to open the file " + std::string(file_path));

   // Get file size
   std::fseek(fp, 0L, SEEK_END);
   int64_t file_size = std::ftell(fp);
   std::rewind(fp); // back to beginning of file

   char *data = (char *)calloc(file_size, 1);

   // here assume that fread does not append extra null-termination (shoudl not)
   size_t n_char = std::fread((void *)data, 1, file_size, fp);
   if (static_cast<int64_t>(n_char) != file_size)
      logger(Logger::ERROR, "[compute_file_checksum]: Read" + std::format("{:d}", n_char)
                                + " character where size of file " + std::string(file_path) + " is "
                                + std::format("{:d}", file_size));

   fclose(fp);

   uint64_t computed_checksum = crc64((uint8_t *)data, file_size);
   free(data);
   return computed_checksum;
}

} // namespace Honeycomb