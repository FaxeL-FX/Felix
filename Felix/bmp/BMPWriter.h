#pragma once
#include <vector>
#include <string>

class BMPWriter
{
public:
    static void write_image(std::vector<std::vector<float>> rgb, int height, int width, char* file_name);
    static int get_index(int heigth, int width, int x, int y);
    static float max_float;
};