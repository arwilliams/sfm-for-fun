#include "image/image.hh"

namespace sfm {
namespace image {

Image::Image(int width, int height, int channels, std::vector<unsigned char> buffer)
    : width_(width), height_(height), channels_(channels),
      buffer_(std::move(buffer)) {}

int Image::width() const { return width_; }
int Image::height() const { return height_; }
int Image::channels() const { return channels_; }

}
}
