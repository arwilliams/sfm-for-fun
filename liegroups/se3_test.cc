#include "liegroups/se3.hh"

#include "util/test_helpers.hh"

#include <cassert>
#include <iostream>

namespace sfm {
namespace liegroups {
namespace {

void assert_near(const SE3 &x, const SE3 &y, const double tolerance) {
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            util::assert_near(x.as_matrix()(i, j), y.as_matrix()(i, j), tolerance);
        }
    }
}

void exp_of_zero_is_identity() {
    const SE3 should_be_identity = SE3::exp(SE3::AlgebraVector::zero());
    constexpr double TOL = 1e-20;
    assert_near(should_be_identity, SE3(), TOL);
}

void log_of_identity_is_zero() {
    const SE3::AlgebraVector should_be_zero = SE3().log();
    for (int i = 0; i < SE3::DOF; ++i) {
        assert(should_be_zero(i) == 0.0);
    }
}

}
}
}

int main() {
    std::cout << "exp_of_zero_is_identity..." << std::endl;
    sfm::liegroups::exp_of_zero_is_identity();

    std::cout << "log_of_identity_is_zero..." << std::endl;
    sfm::liegroups::log_of_identity_is_zero();

    std::cout << "All tests passed!" << std::endl;
    return 0;
}
