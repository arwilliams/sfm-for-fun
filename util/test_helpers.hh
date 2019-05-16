#include <cassert>
#include <cmath>

#include "linalg/matrix.hh"

namespace sfm {
namespace util {

void assert_near(const double a, const double b, const double tolerance) {
    assert(std::abs(a - b) < tolerance);
}

template <class G>
typename G::DifferentialMapping exp_diff_numerical(
        const typename G::DifferentialType &v,
        const double step_size = 1e-5) {
    typename G::DifferentialMapping diff;
    for (int i = 0; i < G::DOF; ++i) {
        const typename G::DifferentialType step =
            step_size * G::DifferentialType::unit_vec(i);
        diff.col(i) = (G::exp(v + step) * (G::exp(v)).inverse()).log() / step_size;
    }
    return diff;
}

}
}
