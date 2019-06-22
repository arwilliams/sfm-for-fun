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

    const SO3::AlgebraTransformation skew = SO3::adjoint(w);

    const linalg::Matrix3d V =
        linalg::Matrix3d::identity() +
        coeffs.b * skew +
        coeffs.c * skew * skew;

    const linalg::Vector3d u = SE3::translational_component(wu);
    SO3::AlgebraTransformation d_exp_w;
    const SE3 x = SE3(SO3::exp(w, d_result_by_input ? &d_exp_w : nullptr), V * u);
    if (d_result_by_input) {
        constexpr double SMALL_THETA = 1e-8;
        const double theta_sq = theta * theta;
        const double K = (theta < SMALL_THETA) ?
            -1./12. + theta_sq / 180.- theta_sq * theta_sq / 6720.
            : (coeffs.a - 2 * coeffs.b) / theta_sq;
        const double L = (theta < SMALL_THETA) ?
            -1/60. + theta_sq / 1260. + theta_sq * theta_sq / 60480.
            : (coeffs.b - 3 * coeffs.c) / theta_sq;
        const linalg::Matrix3d W =
            (coeffs.c - coeffs.b) * linalg::Matrix3d::identity() +
            K * skew +
            L * w * w.transpose();
        d_result_by_input->block<3, 3>(0, 0) = d_exp_w;
        d_result_by_input->block<3, 3>(3, 0) = 
            coeffs.b * SO3::adjoint(u) +
            coeffs.c * (w * u.transpose() + u * w.transpose()) +
            w.dot(u) * W;
        d_result_by_input->block<3, 3>(0, 3) = linalg::Matrix3d::zero();
        d_result_by_input->block<3, 3>(3, 3) = d_exp_w;
    }

    return x;
}

SE3::AlgebraVector SE3::log(AlgebraTransformation *const d_result_by_self) const {

    SO3::AlgebraTransformation d_log_by_w;
    const SO3::AlgebraVector w =
        d_result_by_self ? rotation_.log(&d_log_by_w) : rotation_.log();

    const SO3::AlgebraTransformation skew = SO3::adjoint(w);
    const double theta = w.norm();

    const ExpCoefficients coeffs = compute_exp_coefficients(theta);

    const linalg::Matrix3d Vinv =
        linalg::Matrix3d::identity() -
        0.5 * skew +
        coeffs.e * skew * skew;

    SE3::AlgebraVector wu;
    SE3::rotational_component(wu) = w;
    SE3::translational_component(wu) = Vinv * translation_;

    if (d_result_by_self) {
        constexpr double SMALL_THETA = 1e-8;
        const double theta_sq = theta * theta;
        const double K = (theta < SMALL_THETA) ?
            -1./12. + theta_sq / 180.- theta_sq * theta_sq / 6720.
            : (coeffs.a - 2 * coeffs.b) / theta_sq;
        const double L = (theta < SMALL_THETA) ?
            -1/60. + theta_sq / 1260. + theta_sq * theta_sq / 60480.
            : (coeffs.b - 3 * coeffs.c) / theta_sq;
        const linalg::Matrix3d W =
            (coeffs.c - coeffs.b) * linalg::Matrix3d::identity() +
            K * skew +
            L * w * w.transpose();

        const linalg::Vector3d u = SE3::translational_component(wu);
        const linalg::Matrix3d B =
            coeffs.b * SO3::adjoint(u) +
            coeffs.c * (w * u.transpose() + u * w.transpose()) +
            w.dot(u) * W;

        d_result_by_self->block<3, 3>(0, 0) = d_log_by_w;
        d_result_by_self->block<3, 3>(3, 0) = -d_log_by_w * B * d_log_by_w; 
        d_result_by_self->block<3, 3>(0, 3) = linalg::Matrix3d::zero();
        d_result_by_self->block<3, 3>(3, 3) = d_log_by_w;
    }

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
