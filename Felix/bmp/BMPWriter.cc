#include "BMPWriter.h"
#include <Windows.h>
#include <vector>
#include <string>

float BMPWriter::max_float = 1;

void BMPWriter::write_image(std::vector<std::vector<float>> rgb, int height, int width, char* file_name)
{
    HANDLE file = CreateFileA
    (
            file_name,
            GENERIC_WRITE,
            FILE_SHARE_WRITE,
            nullptr,
            CREATE_ALWAYS,
            FILE_ATTRIBUTE_NORMAL,
            nullptr
    );

    std::string file_cnt;

    file_cnt.append("P3\n");
    file_cnt.append(std::to_string(width) + " " + std::to_string(height) + "\n");
    file_cnt.append("255\n");

    for (int y = height - 1; y >= 0; y--)
        for (int x = 0; x < width; x++)
        {
            int offset = (y * width + x);
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


    DWORD written = 0;
    char* ptr = (char*)file_cnt.c_str();
    WriteFile(file,ptr,file_cnt.size(),&written,NULL);
    CloseHandle(file);
}
int BMPWriter::get_index(int heigth, int width, int x, int y) {
    return y * width + x;
}