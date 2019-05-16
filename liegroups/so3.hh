#pragma once

#include "linalg/matrix.hh"

namespace sfm {
namespace liegroups {

class SO3 {
 public:
    //
    // Element of the Lie algebra so(3),
    // (i.e. an axis-angle vector).
    //
    using DifferentialType = linalg::Vector<double, DOF>;

    //
    // Matrix representing a linear map from so(3) to so(3)
    //
    using DifferentialMapping = linalg::Matrix<double, DOF, DOF>;

    //
    // Defaults to the identity.
    //
    SO3();

    //
    // Does not project the input matrix.
    //
    explicit SO3(const linalg::Matrix3d &rot);

    //
    // SO(3) has dimension 3
    //
    static constexpr int DOF = 3;

    //
    // Group product.
    //
    SO3 operator*(const SO3 &other) const;

    //
    // Group inverse.
    //
    SO3 inverse() const;

    //
    // Group action: left matrix multiplication.
    //
    linalg::Vector3d operator*(const linalg::Vector3d &x) const;

    //
    // Exponential map and its right-invariant differential.
    //
    static SO3 exp(const DifferentialType &w,
                   DifferentialMapping *const d_result_by_input = nullptr);

    //
    // Inverse of exponential and its right-invariant differential.
    //
    DifferentialType log(DifferentialMapping *const d_result_by_self = nullptr) const;

    const linalg::Matrix3d &as_matrix() const { return rot_; }

    //
    // For convenience
    //
    static linalg::Matrix3d skew_matrix(const DifferentialType &w); 

    static DifferentialType axis_angle_from_skew(const linalg::Matrix3d &skew);
 private:
    linalg::Matrix3d rot_;
};

}  // namespace liegroups
}  // namespace sfm
