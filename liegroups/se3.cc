#include "liegroups/se3.hh"

#include "liegroups/exp_coefficients.hh"

namespace sfm {
namespace liegroups {

SE3::SE3() : rotation_(SO3()), translation_(linalg::Vector3d::zero()) {}

SE3::SE3(const SO3 &rotation, const linalg::Vector3d &translation)
    : rotation_(rotation),
      translation_(translation) {}

SE3 SE3::rectified() const {
    return SE3(rotation_.rectified(), translation_);
}

SE3 SE3::operator*(const SE3 &other) const {
    return SE3(rotation_ * other.rotation_,
               rotation_ * other.translation_ + translation_);
}

SE3 SE3::inverse() const {
    return SE3(rotation_.inverse(), -(rotation_.inverse() * translation_));
}

linalg::Vector3d SE3::operator*(const linalg::Vector3d &x) const {
    return rotation_ * x + translation_;
}

SE3 SE3::exp(const AlgebraVector &wu, AlgebraTransformation *const d_result_by_input) {
    const SO3::AlgebraVector w = SE3::rotational_component(wu);
    const double theta = w.norm();
    const ExpCoefficients coeffs = compute_exp_coefficients(theta);

    const linalg::Matrix3d skew = SO3::skew_matrix(w);

    const linalg::Matrix3d V =
        linalg::Matrix3d::identity() +
        coeffs.b * skew +
        coeffs.c * skew * skew;

    const linalg::Vector3d u = SE3::translational_component(wu);
    return SE3(SO3::exp(w), V * u);
}

SE3::AlgebraVector SE3::log(AlgebraTransformation *const d_result_by_self) const {

    SO3::AlgebraTransformation d_log_by_w;
    const SO3::AlgebraVector w =
        d_result_by_self ? rotation_.log(&d_log_by_w) : rotation_.log();

    const linalg::Matrix3d skew = SO3::skew_matrix(w);
    const double theta = w.norm();

    const ExpCoefficients coeffs = compute_exp_coefficients(theta);

    const linalg::Matrix3d Vinv =
        linalg::Matrix3d::identity() -
        0.5 * skew +
        coeffs.e * skew * skew;

    SE3::AlgebraVector wu;
    SE3::rotational_component(wu) = w;
    SE3::translational_component(wu) = Vinv * translation_;

    return wu;
}

linalg::Matrix4d SE3::as_matrix() const {
    linalg::Matrix4d mat = linalg::Matrix4d::identity();
    mat.block<3, 3>(0, 0) = rotation_.as_matrix();
    mat.block<3, 1>(0, 3) = translation_;
    return mat;
}

linalg::Vector3d SE3::apply_action(const linalg::Vector3d &x,
                                   linalg::Matrix<double, 3, DOF> &diff) const {
    const linalg::Vector3d transformed = (*this) * x;
    diff.block<3, 3>(0, 0) = SO3::adjoint(-transformed);
    diff.block<3, 3>(0, 3) = linalg::Matrix3d::identity();
    return transformed;
}

}
}
