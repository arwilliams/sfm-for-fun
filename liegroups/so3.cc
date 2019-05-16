#include "liegroups/so3.hh"

#include <cassert>
#include <cmath>

#include "math/clamp.hh"

namespace sfm {
namespace liegroups {
namespace {

constexpr double SMALL_X = 1e-8;

double exp_coeff_a(const double x) {
	if (x < SMALL_X) {
		const double x_sq = x * x;
		return 1. - x_sq / 6. + x_sq * x_sq / 120.;
	}
	return std::sin(x) / x;
}

double exp_coeff_b(const double x) {
	if (x < SMALL_X) {
		const double x_sq = x * x;
		return 0.5 - x_sq / 24. + x_sq * x_sq / 720.;
	}
	return (1 - std::cos(x)) / (x * x);
}

double exp_coeff_c(const double x) {
    if (x < SMALL_X) {
        const double x_sq = x * x;
        return 1. / 6. - x_sq / 120. + x_sq * x_sq / 5040.;
    }
    return (1. - exp_coeff_a(x)) / (x * x);
}

//void project_onto_SO3(linalg::Matrix3d &mat) {
//    // TODO
//    return;
//}
//
}

SO3::SO3() : rot_(linalg::Matrix3d::identity()) {}

SO3::SO3(const linalg::Matrix3d &rot)
    : rot_(rot) {}

SO3 SO3::operator*(const SO3 &other) const {
    return SO3(rot_ * other.rot_);
}

SO3 SO3::exp(const DifferentialType &w,
             DifferentialMapping *const d_result_by_input) {
    const double theta = w.norm();
    const double a = exp_coeff_a(theta);
    const double b = exp_coeff_b(theta);

    const linalg::Matrix3d skew = SO3::skew_matrix(w);

    if (d_result_by_input) {
        const double c = exp_coeff_c(theta);
        *d_result_by_input = a * linalg::Matrix3d::identity() +
                             b * skew +
                             c * w * w.transpose();
    }

    return SO3(linalg::Matrix3d::identity() + a * skew + b * skew * skew);
}

SO3::DifferentialType SO3::log(SO3::DifferentialMapping *const d_result_by_input) const {

    // Extract the skew-symmetric part w --- 
    // to get the normalized axis we need to divde w by sin(theta),
    // then the correct norm is theta.
    SO3::DifferentialType w =
        SO3::axis_angle_from_skew(0.5 * (rot_ - rot_.transpose()));

    const double trace = rot_.trace();
    const double cos_theta = 0.5 * (trace - 1.0);
    const double sin_theta_sq = w.squared_norm();

    if (cos_theta > 1. - 1e-8) {
        //
        // Small angle --- sin(theta) / theta close to 1, so
        // there is enough information in the anti-symmetric part.
        //
        // We can evaluate the Taylor series for arcsin(x) / x
        // for x = sin(theta), yielding theta / sin(theta).
        // uses only even powers so can compute using sin_theta_sq,
        // without taking square root.
        //
        const double inv_sin_theta =
            1.0 + sin_theta_sq / 6.0
                + sin_theta_sq * sin_theta_sq * (3.0 / 40.0)
                + sin_theta_sq * sin_theta_sq * sin_theta_sq * (5.0 / 112.0);

        return w * inv_sin_theta;
    }

    if (cos_theta > -1. + 1e-8) {
        // Reasonable angle.
        const double theta = std::acos(cos_theta);
        const double sin_theta = std::sqrt(std::max(sin_theta_sq, 0.));
        return (theta / sin_theta) * w;
    }

    //
    // Angle near pi --- sin(theta) / theta is close to zero, so
    // we can't use the symmetric part. instead we must factor
    // antisymmetric part, making an arbitrary choice of sign in the
    // process.
    //
    const double sin_theta = std::sqrt(std::max(sin_theta_sq, 0.));
    const double theta = M_PI - std::asin(math::clamp(sin_theta, -1., 1.));

    const linalg::Matrix3d B = 0.5 * (rot_ + linalg::Matrix3d::identity());
    const double a_sq = B(0, 0);
    const double abs_a = std::sqrt(a_sq);

    const double b_sq = B(1, 1);
    const double abs_b = std::sqrt(b_sq);

    const double c_sq = B(2, 2);
    const double abs_c = std::sqrt(c_sq);

    const double ab = B(0, 1);
    short sgn_a, sgn_b, sgn_c;
    sgn_a = 1;
    if (ab < 0) {
        sgn_b = -1;
    } else {
        sgn_b = 1;
    }

    const double ac = B(0, 2);
    if (ac < 0) {
        sgn_c = -1;
    } else {
        sgn_c = 1;
    }

    return theta * DifferentialType(abs_a * sgn_a, abs_b * sgn_b, abs_c * sgn_c);
}

SO3 SO3::inverse() const {
    return SO3(rot_.transpose());
}

linalg::Matrix3d SO3::skew_matrix(const DifferentialType &w) {
    linalg::Matrix3d wx = linalg::Matrix3d::zero();
    wx(1, 0) = w(2, 0);
    wx(0, 1) = -w(2, 0);

    wx(2, 0) = -w(1, 0);
    wx(0, 2) = w(1, 0);

    wx(2, 1) = w(0, 0);
    wx(1, 2) = -w(0, 0);

    return wx;
}

SO3::DifferentialType SO3::axis_angle_from_skew(const linalg::Matrix3d &skew) {
    return DifferentialType(skew(2, 1), skew(0, 2), skew(1, 0));
}

}  // namespace liegroups
}  // namespace sfm
