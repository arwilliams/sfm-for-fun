#pragma once

#include "linalg/matrix.hh"

namespace sfm {
namespace geometry {

template <int DIM>
inline linalg::Vector<double, DIM> project(
    const linalg::Vector<double, DIM> &to_project,
    const linalg::Vector<double, DIM> &onto) {
    return to_project.dot(onto) / onto.dot(onto) * onto;
}

}
}
