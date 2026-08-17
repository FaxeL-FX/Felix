#include "BMPWriter.h"
#include <vector>
#include <string>

float BMPWriter::max_float = 1;

std::string BMPWriter::make_bmp_image(const std::vector<std::vector<float>> &rgb, int size)
{
    std::string file_cnt;

    file_cnt.append("P3\n");
    file_cnt.append(std::to_string(size) + " " + std::to_string(size) + "\n");
    file_cnt.append("255\n");

    for (int y = size - 1; y >= 0; y--)
        for (int x = 0; x < size; x++)
        {
            int offset = (y * size + x);
            if ((rgb[offset]).size() < 3) {
                file_cnt.append("0 0 0\n");
                continue;
            }
            int val_red = (int)((rgb[offset][0]) * (float)255);
            int val_green = (int)((rgb[offset][1]) * (float)255);
            int val_blue = (int)((rgb[offset][2]) * (float)255);
            if (val_red < 0)   val_red = 0;
            if (val_green < 0) val_green = 0;
            if (val_blue < 0)  val_blue = 0;
            if (255 < val_red)   val_red = 255;
            if (255 < val_green) val_green = 255;
            if (255 < val_blue)  val_blue = 255;
            file_cnt.append(std::to_string(val_red) + " " + std::to_string(val_green) + " " + std::to_string(val_blue) + "\n");
        }

    /*for (int y = height - 1; y >= 0; y--)
        for (int x = 0; x < width; x++)
        {
            int offset = (y * width + x);
            int val_red = (int)((red[offset]) * ((float)255 / max_float));
            int val_green = (int)((green[offset]) * ((float)255 / max_float));
            int val_blue = (int)((blue[offset]) * ((float)255 / max_float));
            if (val_red   < 0)   val_red = 0;
            if (val_green < 0) val_green = 0;
            if (val_blue  < 0)  val_blue = 0;
            file_cnt.append(std::to_string(val_red) + " " + std::to_string(val_green) + " " + std::to_string(val_blue) + "\n");
        }*/

    return file_cnt;
}
int BMPWriter::get_index(std::vector<std::vector<float>>* rgb, int x, int y) {
    return y * rgb->size() + x;
}