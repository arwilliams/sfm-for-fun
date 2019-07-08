#pragma once

#include "image/image.hh"

#include <string>

namespace sfm {
namespace image {

Image load_image_from_png(const std::string &filename);

}
}
