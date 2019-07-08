#include "image/load_image_from_png.hh"

#include <cassert>
#include <stdexcept>

#include "stdio.h"
#include "setjmp.h"

#include "png.h"

namespace sfm {
namespace image {

//
// http://www.libpng.org/pub/png/book/chapter13.html#png.ch13.div.1
//
Image load_image_from_png(const std::string &filename) {
    constexpr char MODE[] = "r";
    FILE *fin = fopen(filename.c_str(), MODE);
    if (!fin) {
        fclose(fin);
        throw std::runtime_error("Failed to open file for reading: " + filename);
    }

    unsigned char sig[8];
    fread(sig, 1, 8, fin);
    if (!png_check_sig(sig, 8)) {
        fclose(fin);
        throw std::runtime_error("png_check_sig failed: " + filename);
    }

    png_structp png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING,
                                                 NULL,
                                                 NULL,
                                                 NULL);
    if (!png_ptr) {
        fclose(fin);
        throw std::runtime_error("png_create_read_struct failed: " + filename);
    }

    png_infop info_ptr = png_create_info_struct(png_ptr);
    if (!info_ptr) {
        png_destroy_read_struct(&png_ptr, NULL, NULL);
        fclose(fin);
        throw std::runtime_error("png_create_info_struct failed: " + filename);
    }

    // https://github.com/rubis-lab/CPSim_linux/issues/1
    if (setjmp(png_jmpbuf(png_ptr))) {
        png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
        fclose(fin);
        throw std::runtime_error("setjmp on png ptr jmpbuf failed: " + filename);
    }

    png_init_io(png_ptr, fin);
    png_set_sig_bytes(png_ptr, 8);
    png_read_info(png_ptr, info_ptr);

    png_uint_32 width, height;
    int bit_depth, color_type;
    if (!png_get_IHDR(png_ptr, info_ptr, &width, &height, &bit_depth, &color_type,
                      NULL, NULL, NULL)) {
        png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
        fclose(fin);
        throw std::runtime_error("png_get_IHDR failed: " + filename);
    }

    if (setjmp(png_jmpbuf(png_ptr))) {
        png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
        fclose(fin);
        throw std::runtime_error("setjmp on png ptr jmpbuf failed: " + filename);
    }

    const png_uint_32 row_bytes = png_get_rowbytes(png_ptr, info_ptr);
    const int channels = png_get_channels(png_ptr, info_ptr);

    assert(channels * width == row_bytes);

    std::vector<unsigned char> image_data(row_bytes * height);

    png_bytep row_pointers[height];
    for (uint32_t i = 0; i < height; ++i) {
        row_pointers[i] = image_data.data() + i * row_bytes;
    }

    png_read_image(png_ptr, row_pointers);

    const Image image(width, height, channels, std::move(image_data));

    png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
    fclose(fin);
    return image;
}

}
}
