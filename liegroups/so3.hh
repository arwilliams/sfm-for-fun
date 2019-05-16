#pragma once

#include "linalg/matrix.hh"

namespace sfm {
namespace liegroups {

class SO3 {
 public:
    SO3();

    explicit SO3(const linalg::Matrix3d &rot);

    static constexpr int DOF = 3;

    using DifferentialType = linalg::Vector<double, DOF>;

    using DifferentialMapping = linalg::Matrix<double, DOF, DOF>;

    SO3 operator*(const SO3 &other) const;

    static SO3 exp(const DifferentialType &w,
                   DifferentialMapping *const d_result_by_input = nullptr);

    DifferentialType log(DifferentialMapping *const d_result_by_self = nullptr) const;

    SO3 inverse() const;

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
