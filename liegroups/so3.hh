#pragma once

#include "linalg/matrix.hh"

namespace sfm {
namespace liegroups {

//
// An element of the Lie group SO(3), representing a
// rigid transformation of 3-space preserving the origin.
//
class SO3 {
 public:
    //
    // SO(3) has dimension 3
    //
    static constexpr int DOF = 3;

    //
    // Element of the Lie algebra so(3),
    // (i.e. an axis-angle vector).
    //
    using AlgebraVector = linalg::Vector<double, DOF>;

    //
    // Matrix representing a linear map from so(3) to so(3)
    //
    using AlgebraTransformation = linalg::Matrix<double, DOF, DOF>;

    //
    // Defaults to the identity.
    //
    SO3();

    //
    // Does not project the input matrix.
    //
    explicit SO3(const linalg::Matrix3d &rot);

    //
    // Returns an SO(3) whose underlying matrix is the nearest
    // orthogonal matrix to this one's.
    //
    SO3 rectified() const;

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
    static SO3 exp(const AlgebraVector &w,
                   AlgebraTransformation *const d_result_by_input = nullptr);

    //
    // Inverse of exponential and its right-invariant differential.
    //
    AlgebraVector log(AlgebraTransformation *const d_result_by_self = nullptr) const;

    //
    // Adjoint representation of a group element
    //
    AlgebraTransformation adjoint() const;

    //
    // Adjoint representation of a Lie algebra element
    //
    static AlgebraTransformation adjoint(const AlgebraVector &w);

    const linalg::Matrix3d &as_matrix() const { return rot_; }

    //
    // Applies the group action per operator* and returns the
    // right-invariant differential with respect to the group.
    //
    linalg::Vector3d apply_action(const linalg::Vector3d &x,
                                  linalg::Matrix<double, 3, DOF> &diff) const;

    //
    // For convenience
    //
    static linalg::Matrix3d skew_matrix(const AlgebraVector &w); 

    static AlgebraVector axis_angle_from_skew(const linalg::Matrix3d &skew);
 private:
    linalg::Matrix3d rot_;
};

}  // namespace liegroups
}  // namespace sfm
