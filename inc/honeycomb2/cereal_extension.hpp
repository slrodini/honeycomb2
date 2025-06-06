#ifndef HC2_CERAL_EXTENSION_HPP
#define HC2_CERAL_EXTENSION_HPP

#include <honeycomb2/cereal/archives/portable_binary.hpp>
#include <honeycomb2/cereal/types/vector.hpp>
#include <honeycomb2/cereal/types/utility.hpp>
#include <honeycomb2/cereal/types/string.hpp>
#include <honeycomb2/cereal/access.hpp>
#include <honeycomb2/cereal/types/memory.hpp>
#include <honeycomb2/cereal/details/traits.hpp>

#include <sstream>
#include <fstream>
#include <vector>
#include <honeycomb2/utilities.hpp>
#include <honeycomb2/checksum.hpp>

/**
 * @file cereal_extension.hpp
 * @author Simone Rodini (rodini.simone.luigi@gmail.com)
 * @brief  Functions to interface with cereal and to incorporate checksum into the archives.
 * @version 0.1
 * @date 2025-06-04
 *
 * @copyright Copyright (c) 2025
 *
 */

namespace Honeycomb
{
/**
 *
 * @brief Archives data using cereal, appending the checksum to the file
 *
 * @tparam T        The type that is being serialized.
 * @tparam Archive  The archive for serialization.
 * @param data      The data to be serialized.
 * @param file_name The name of the output file
 */
template <typename T, class Archive>
void SaveChecksumArchive(const T &data, const std::string &file_name)
{
   if constexpr (!std::is_base_of<cereal::detail::OutputArchiveBase, Archive>::value) {
      logger(Logger::ERROR,
             "Trying to use `SaveChecksumArchive` with an archive that is not an output archive.");
   }
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

/**
 * @brief Load data from disk and verify their integrity via the checksum
 *
 * @tparam T         The type of data to be read.
 * @tparam Archive   The archive format to be used.
 * @param file_name  The name of the file to load.
 * @param data       The destination for the data to be read.
 * @return True on successful load, False on faulty load.
 */
template <typename T, class Archive>
bool LoadAndVerify(const std::string &file_name, T &data)
{
   if constexpr (!std::is_base_of<cereal::detail::InputArchiveBase, Archive>::value) {
      logger(Logger::ERROR,
             "Trying to use `LoadAndVerify` with an archive that is not an input archive.");
   }
   // std::ifstream file(file_name, std::ios::binary | std::ios::ate);
   std::FILE *file = std::fopen(file_name.c_str(), "rb");

   if (nullptr == file) {
      logger(Logger::WARNING, "[LoadAndVerify]: Failed to open file: " + file_name);
      return false;
   }

   std::fseek(file, 0L, SEEK_END);
   uint64_t file_size = std::ftell(file);
   std::rewind(file);

   if (file_size < sizeof(uint64_t))
      logger(Logger::ERROR, "[LoadAndVerify]: File: " + file_name + "too short");

   std::vector<char> buffer(file_size);
   size_t n_char = std::fread(buffer.data(), 1, file_size, file);
   std::fclose(file);

   if (n_char != file_size)
      logger(Logger::ERROR, "[LoadAndVerify]: Read" + std::format("{:d}", n_char)
                                + " character where size of file " + file_name + " is "
                                + std::format("{:d}", file_size));
   // Compute checksum of data
   uint64_t computed_checksum = crc64((const uint8_t *)buffer.data(), file_size);

   // Verify checksum
   if (computed_checksum)
      logger(Logger::ERROR, "[LoadAndVerify]: File: " + file_name + "has incorrect checksum! Got: "
                                + std::format("{:d}", computed_checksum));

   std::stringstream ss;
   ss.write(buffer.data(), file_size - sizeof(uint64_t));
   Archive ar(ss);
   ar(data);
   return true;
}

/**
 * @brief Same as \ref Honeycomb::SaveChecksumArchive, but without checksum
 *
 * @tparam T        The type that is being serialized.
 * @tparam Archive  The archive for serialization.
 * @param data      The data to be serialized.
 * @param file_name The name of the output file
 *
 */
template <typename T, class Archive>
void SaveArchive(const T &data, const std::string &file_name)
{
   if constexpr (!std::is_base_of<cereal::detail::OutputArchiveBase, Archive>::value) {
      logger(Logger::ERROR,
             "Trying to use `SaveArchive` with an archive that is not an output archive.");
   }
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

/**
 * @brief Same as \ref Honeycomb::LoadAndVerify but without checksum
 *
 * @tparam T         The type of data to be read.
 * @tparam Archive   The archive format to be used.
 * @param file_name  The name of the file to load.
 * @param data       The destination for the data to be read.
 * @return True on successful load, False on faulty load.
 *
 */
template <typename T, class Archive>
bool LoadArchive(const std::string &file_name, T &data)
{
   if constexpr (!std::is_base_of<cereal::detail::InputArchiveBase, Archive>::value) {
      logger(Logger::ERROR,
             "Trying to use `LoadArchive` with an archive that is not an input archive.");
   }
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
      logger(Logger::ERROR, " Read" + std::format("{:d}", n_char) + " character where size of file "
                                + file_name + " is " + std::format("{:d}", file_size));

   std::stringstream ss;
   ss.write(buffer.data(), file_size);
   Archive ar(ss);
   ar(data);
   return true;
}

} // namespace Honeycomb

#endif // HC2_CERAL_EXTENSION_HPP
