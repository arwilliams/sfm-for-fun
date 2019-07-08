#include "image/load_image_from_png.hh"

#include <cassert>
#include <iostream>

struct ImageInfo {
    int width;
    int height;
    int channels;
};

ImageInfo get_info_for_known_image(const std::string &filepath) {
    if (filepath.find("pngkey.com-cat-silhouette-png-1527338.png") != std::string::npos) {
        return ImageInfo{
            .width = 1920,
            .height = 1155,
            .channels = 4,
        };
    }
    std::cerr << "Unknown image: " << filepath << std::endl;
    assert(false);
}

int main(int argc, char *argv[]) {
    if (argc != 2) {
        std::cerr << "Image path required." << std::endl;
        return 1;
    }

    const char *image_path = argv[1];

    const sfm::image::Image image = sfm::image::load_image_from_png(std::string(image_path));

    const ImageInfo info = get_info_for_known_image(std::string(image_path));
    assert(image.width() == info.width);
    assert(image.height() == info.height);
    assert(image.channels() == info.channels);

    std::cout << "All tests passed!" << std::endl;
    return 0;
}
