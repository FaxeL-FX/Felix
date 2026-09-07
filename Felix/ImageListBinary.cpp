#define STB_IMAGE_WRITE_IMPLEMENTATION

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iterator>
#include <vector>

#include "bmp/gif.h"
#include "bmp/stb_image_write.h"

#include "Felix.h"
#include "ImageListBinary.hh"

static void write_gif_frame(
    GifWriter *writer,
    const std::vector<std::vector<float>> &img,
    int resolution,
    uint32_t delay)
{
    std::vector<uint8_t> image;
    image.resize(resolution * resolution * 4);
    for (int yy = 0; yy < resolution; ++yy) {
        for (int xx = 0; xx < resolution; ++xx) {
            size_t offset_image = ((resolution - yy - 1) * resolution + xx);
            size_t offset = (yy * resolution + xx);
            auto &r = image[offset_image * 4 + 0];
            auto &g = image[offset_image * 4 + 1];
            auto &b = image[offset_image * 4 + 2];
            auto &a = image[offset_image * 4 + 3];
            r = img[offset][0] * 255;
            g = img[offset][1] * 255;
            b = img[offset][2] * 255;
            a = 255;
        }
    }
    GifWriteFrame(writer, image.data(), resolution, resolution, delay);
}

ImageListBinary::ImageListBinary(const Felix::ImageList &images)
{
    if (images.empty()) {
        return;
    }

    std::vector<uint8_t> &output = *static_cast<std::vector<uint8_t> *>(this);

    // png
    if (images.size() == 1) {
        auto &image = images[0];

        std::vector<uint8_t> rgb;
        rgb.resize(image.resolution * image.resolution * 3);
        for (size_t y = 0; y != image.resolution; ++y) {
            for (size_t x = 0; x != image.resolution; ++x) {
                size_t ry = image.resolution - y - 1;
                size_t offset = (y * image.resolution + x);
                auto &pixel = image[offset];
                if (pixel.size() != 3) {
                    continue;
                }
                size_t roffset = (ry * image.resolution + x);
                size_t rgb_offset = roffset * 3;
                rgb[rgb_offset + 0] = pixel[0] * 255;
                rgb[rgb_offset + 1] = pixel[1] * 255;
                rgb[rgb_offset + 2] = pixel[2] * 255;
            }
        }

        int size = 0;
        auto png = stbi_write_png_to_mem(
            rgb.data(),
            image.resolution * 3,
            image.resolution,
            image.resolution,
            3,
            &size);
        try {
            std::copy(png, png + size, std::back_inserter(output));
        } catch (...) {
            free(png);
            throw;
        }
        free(png);
        file_name = "graph.png";
        return;
    }

    // gif
    {
        auto resolution = images[0].resolution;
        std::vector<uint8_t> data;
        GifWriterCallback cb;
        cb.udata = &output;
        cb.write = [](void *udata, void *data, size_t size) {
            std::vector<uint8_t> &vec =
                *static_cast<std::vector<uint8_t> *>(udata);
            uint8_t *start = (uint8_t *)data;
            uint8_t *end = start + size;
            std::copy(start, end, std::back_inserter(vec));
        };

        GifWriter writer = {};
        GifBegin(&writer, cb, resolution, resolution, images.gif_delay);
        for (auto &image : images) {
            write_gif_frame(&writer, image, resolution, images.gif_delay);
        }
        GifEnd(&writer);
        file_name = "graph.gif";
        return;
    }
}
