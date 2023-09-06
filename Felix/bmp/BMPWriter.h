#pragma once
#include <vector>
#include <string>

class BMPWriter
{
public:
    static void write_image(std::vector<std::vector<float>> rgb, int size, char* file_name);
    static int get_index(std::vector<std::vector<float>>* rgb, int x, int y);
    static float max_float;
};