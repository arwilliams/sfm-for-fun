#pragma once

#include <vector>

namespace sfm {
namespace image {

class Image {
 public:
    Image() = default;
    Image(int width,
          int height,
          int channels,
          std::vector<unsigned char> buffer);

    int width() const;
    int height() const;
    int channels() const;

 private:
    int width_;
    int height_;
    int channels_;

    std::vector<unsigned char> buffer_;
};

}
}
