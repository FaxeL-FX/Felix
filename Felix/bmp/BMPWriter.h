#pragma once
#include <vector>
#include <string>

class BMPWriter
{
public:
    static std::string make_bmp_image(const std::vector<std::vector<float>> &rgb, int size);
    static int get_index(std::vector<std::vector<float>>* rgb, int x, int y);
    static float max_float;
};