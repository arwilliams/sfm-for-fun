#pragma once

#include <cassert>
#include <cmath>

#include "linalg/matrix.hh"

namespace sfm {
namespace util {

void assert_near(const double a, const double b, const double tolerance) {
    assert(std::abs(a - b) < tolerance);
}

template <class G>
typename G::AlgebraTransformation exp_diff_numerical(
        const typename G::AlgebraVector &v,
        const double step_size = 1e-5) {
    typename G::AlgebraTransformation diff;
    for (int i = 0; i < G::DOF; ++i) {
        const typename G::AlgebraVector step =
            step_size * G::AlgebraVector::unit_vec(i);
        diff.col(i) = (G::exp(v + step) * (G::exp(v)).inverse()).log() / step_size;
    }
    return diff;
}

template <class G>
typename G::AlgebraTransformation log_diff_numerical(
    const G &x,
    const double step_size = 1e-5) {
    typename G::AlgebraTransformation diff;
    for (int i = 0; i < G::DOF; ++i) {
        const typename G::AlgebraVector step =
            step_size * G::AlgebraVector::unit_vec(i);
        diff.col(i) = ((G::exp(step) * x).log() - x.log()) / step_size;
    }
    return diff;
}

}
}
