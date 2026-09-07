#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "Felix.h"

struct ImageListBinary : std::vector<uint8_t>
{
    ImageListBinary(const Felix::ImageList &images);
    std::string file_name;
};
