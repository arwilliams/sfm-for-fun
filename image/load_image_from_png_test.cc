#include "image/load_image_from_png.hh"

#include <iostream>

int main(int argc, char *argv[]) {
    if (argc != 2) {
        std::cerr << "Image path required." << std::endl;
        return 1;
    }

    const char *image_path = argv[1];

    sfm::image::load_image_from_png(std::string(image_path));

    std::cout << "All tests passed!" << std::endl;
    return 0;
}
