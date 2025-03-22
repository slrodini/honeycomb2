#include <cereal/archives/portable_binary.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/access.hpp>
#include <cereal/types/memory.hpp>

#include <sstream>
#include <fstream>
#include <vector>
#include <honeycomb2/utilities.hpp>
#include <honeycomb2/checksum.hpp>

namespace Honeycomb
{

template <typename T, class Archive>
void SaveChecksumArchive(const T &data, const std::string &file_name)
{
   std::ostringstream buffer;
   uint64_t checksum;
   {
      Archive ar(buffer);
      ar(data);
   }

   std::string serialized = buffer.str();
   checksum               = crc64((const uint8_t *)serialized.c_str(), serialized.size());

   std::ofstream file(file_name, std::ios::binary);

   file.write(serialized.c_str(), serialized.size());

   // Correctly append the checksum so that the checksum of the file is 0
   for (int i = 0; i < 8; ++i) {
      uint8_t c = (checksum >> (56 - 8 * i)) & 0xFF;
      file.write(reinterpret_cast<const char *>(&c), sizeof(uint8_t));
   }

   file.close();
   checksum = 0;
   buffer.str("");
   buffer.clear();
}

template <typename T, class Archive>
void SaveArchive(const T &data, const std::string &file_name)
{
   std::ostringstream buffer;
   {
      Archive ar(buffer);
      ar(data);
   }
   std::string serialized = buffer.str();
   std::ofstream file(file_name, std::ios::binary);
   file.write(serialized.c_str(), serialized.size());
   file.close();
}

template <typename T, class Archive>
bool LoadAndVerify(const std::string &file_name, T &data)
{
   // std::ifstream file(file_name, std::ios::binary | std::ios::ate);
   std::FILE *file = std::fopen(file_name.c_str(), "rb");

   if (nullptr == file) logger(Logger::ERROR, "[LoadAndVerify]: Failed to open file: " + file_name);

   std::fseek(file, 0L, SEEK_END);
   uint64_t file_size = std::ftell(file);
   std::rewind(file);

   if (file_size < sizeof(uint64_t)) logger(Logger::ERROR, "[LoadAndVerify]: File: " + file_name + "too short");

   std::vector<char> buffer(file_size);
   size_t n_char = std::fread(buffer.data(), 1, file_size, file);
   std::fclose(file);

   if (n_char != file_size)
      logger(Logger::ERROR, "[LoadAndVerify]: Read" + std::format("{:d}", n_char) + " character where size of file " +
                                file_name + " is " + std::format("{:d}", file_size));
   // Compute checksum of data
   uint64_t computed_checksum = crc64((const uint8_t *)buffer.data(), file_size);

   // Verify checksum
   if (computed_checksum)
      logger(Logger::ERROR, "[LoadAndVerify]: File: " + file_name +
                                "has incorrect checksum! Got: " + std::format("{:d}", computed_checksum));

   std::stringstream ss;
   ss.write(buffer.data(), file_size - sizeof(uint64_t));
   Archive ar(ss);
   ar(data);
   return true;
}

template <typename T, class Archive>
bool LoadArchive(const std::string &file_name, T &data)
{
   // std::ifstream file(file_name, std::ios::binary | std::ios::ate);
   std::FILE *file = std::fopen(file_name.c_str(), "rb");

   if (nullptr == file) logger(Logger::ERROR, "Failed to open file: " + file_name);

   std::fseek(file, 0L, SEEK_END);
   uint64_t file_size = std::ftell(file);
   std::rewind(file);

   if (file_size < sizeof(uint64_t)) logger(Logger::ERROR, "File: " + file_name + "too short");

   std::vector<char> buffer(file_size);
   size_t n_char = std::fread(buffer.data(), 1, file_size, file);
   std::fclose(file);

   if (n_char != file_size)
      logger(Logger::ERROR, " Read" + std::format("{:d}", n_char) + " character where size of file " + file_name +
                                " is " + std::format("{:d}", file_size));

   std::stringstream ss;
   ss.write(buffer.data(), file_size);
   Archive ar(ss);
   ar(data);
   return true;
}

} // namespace Honeycomb