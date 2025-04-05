#pragma once

#include <honeycomb2/default.hpp>
#include <honeycomb2/utilities.hpp>

namespace Honeycomb
{

uint64_t crc64(const uint8_t *data, size_t length);
uint64_t compute_file_checksum(const char *file_path);

} // namespace Honeycomb
