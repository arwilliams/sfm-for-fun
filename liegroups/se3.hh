#pragma once

#include "linalg/matrix.hh"
#include "liegroups/so3.hh"

namespace sfm {
namespace liegroups {

//
// An element of the Lie group SE(3), representing a
// rigid transformation of 3-space.
//
class SE3 {
 public:
    //
    // SE(3) has dimension 6
    //
    static constexpr int DOF = 6;

    //
    // Element of the Lie algebra se(3).
    //
    // First three components are rotational; last three
    // are translational.
    //
    using AlgebraVector = linalg::Vector<double, DOF>;

    //
    // Matrix representing a linear map from se(3) to se(3)
    //
    using AlgebraTransformation = linalg::Matrix<double, DOF, DOF>;

    //
    // Defaults to the identity.
    //
    SE3();

    //
    //
    // Constructs the SE(3) element whose action is x -> Rx + t.
    //
    // Does not project the input rotation.
    //
    SE3(const SO3 &rotation, const linalg::Vector3d &translation);

    //
    // Returns an SE(3) whose underlying rotation matrix is the nearest
    // orthogonal matrix to this one's.
    //
    SE3 rectified() const;

    //
    // Group product.
    //
    SE3 operator*(const SE3 &other) const;

    //
    // Group inverse.
    //
    SE3 inverse() const;

    //
    // Group action: left matrix multiplication on affine points.
    //
    linalg::Vector3d operator*(const linalg::Vector3d &x) const;

    //
    // Exponential map and its right-invariant differential.
    //
    static SE3 exp(const AlgebraVector &wu,
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

    linalg::Matrix4d as_matrix() const;

    //
    // Applies the group action per operator* and returns the
    // right-invariant differential with respect to the group.
    //
    linalg::Vector3d apply_action(const linalg::Vector3d &x,
                                  linalg::Matrix<double, 3, DOF> &diff) const;

    //
    // For convenience.
    //
    static const auto rotational_component(const AlgebraVector &wu) {
        return wu.block<3, 1>(0, 0);
    }

    static auto rotational_component(AlgebraVector &wu) {
        return wu.block<3, 1>(0, 0);
    }

    static const auto translational_component(const AlgebraVector &wu) {
        return wu.block<3, 1>(3, 0);
    }

    static auto translational_component(AlgebraVector &wu) {
        return wu.block<3, 1>(3, 0);
    }
 private:
    SO3 rotation_;
    linalg::Vector3d translation_;
};

}  // namespace liegroups
}  // namespace sfm
