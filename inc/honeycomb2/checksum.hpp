#ifndef HC2_CHECKSUM_HPP
#define HC2_CHECKSUM_HPP

#include <honeycomb2/default.hpp>
#include <honeycomb2/utilities.hpp>

/**
 * @file checksum.hpp
 * @author Simone Rodini (rodini.simone.luigi@gmail.com)
 * @brief  Utilities to compute the checksum of data and files
 * @version 0.1
 * @date 2025-06-04
 *
 * @copyright Copyright (c) 2025
 *
 */

namespace Honeycomb
{

/**
 * @brief Computes the CRC64 checksum on a chunk of data
 *
 * @param data      The chunk of data to be processed.
 * @param length    The length (numer of bytes) of the data.
 * @return uint64_t The CRC64 checksum
 *
 * CRC-64/WE, without final xor
 * Resources: https://www.wolfgang-ehrhardt.de/ ;
 * https://reveng.sourceforge.io/crc-catalogue/all.htm
 */
uint64_t crc64(const uint8_t *data, size_t length);
/**
 * @brief Computes the checksum of a file, loaded from disk
 *
 * @param file_path The path to the file.
 * @return uint64_t The CRC64 checksum
 */
uint64_t compute_file_checksum(const char *file_path);

} // namespace Honeycomb

#endif // HC2_
HC2_CHECKSUM_HPP