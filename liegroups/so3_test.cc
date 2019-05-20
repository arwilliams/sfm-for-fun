#include "liegroups/so3.hh"

#include "util/test_helpers.hh"

#include <cassert>
#include <iostream>

namespace sfm {
namespace liegroups {
namespace {

void assert_near(const SO3 &x, const SO3 &y, const double tolerance) {
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            util::assert_near(x.as_matrix()(i, j), y.as_matrix()(i, j), tolerance);
        }
    }
}

SO3 make_rotation_about_x(const double theta) {
    linalg::Matrix3d rot_x = linalg::Matrix3d::identity();
    rot_x(1, 1) = std::cos(theta); rot_x(1, 2) = -std::sin(theta);
    rot_x(2, 1) = std::sin(theta); rot_x(2, 2) = std::cos(theta);
    return SO3(rot_x);
}

SO3 make_rotation_about_y(const double theta) {
    linalg::Matrix3d rot_y = linalg::Matrix3d::identity();
    rot_y(0, 0) = std::cos(theta); rot_y(0, 2) = std::sin(theta);
    rot_y(2, 0) = -std::sin(theta); rot_y(2, 2) = std::cos(theta);
    return SO3(rot_y);
}

SO3 make_rotation_about_z(const double theta) {
    linalg::Matrix3d rot_z = linalg::Matrix3d::identity();
    rot_z(0, 0) = std::cos(theta); rot_z(0, 1) = -std::sin(theta);
    rot_z(1, 0) = std::sin(theta); rot_z(1, 1) = std::cos(theta);
    return SO3(rot_z);
}

}

void exp_of_zero_is_identity() {
    const linalg::Matrix3d should_be_identity =
        SO3::exp(SO3::AlgebraVector::zero()).as_matrix();
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            if (i == j) {
                assert(should_be_identity(i, j) == 1.);
            } else {
                assert(should_be_identity(i, j) == 0.);
            }
        }
    }
}

void exp_of_rotation_around_x() {
    constexpr double THETA_X = M_PI / 8;
    const SO3::AlgebraVector around_x_axis(THETA_X, 0., 0.);
    const SO3 rot_x = SO3::exp(around_x_axis);

    const SO3 true_rot_x = make_rotation_about_x(THETA_X);
    constexpr double TOL = 1e-10;
    assert_near(rot_x, true_rot_x, TOL);
}

void exp_of_rotation_around_y() {
    constexpr double THETA_Y = -3 * M_PI / 5;
    const SO3::AlgebraVector around_y_axis(0., THETA_Y, 0.);
    const SO3 rot_y = SO3::exp(around_y_axis);
    const SO3 true_rot_y = make_rotation_about_y(THETA_Y);
    constexpr double TOL = 1e-10;
    assert_near(rot_y, true_rot_y, TOL);
}

void exp_of_rotation_around_z() {
    constexpr double THETA_Z = 20 * M_PI / 3;
    const SO3::AlgebraVector around_z_axis(0., 0., THETA_Z);
    const SO3 rot_z = SO3::exp(around_z_axis);
    constexpr double TOL = 1e-10;
    const SO3 true_rot_z = make_rotation_about_z(THETA_Z);
    assert_near(rot_z, true_rot_z, TOL);
}

void exp_of_small_angle_rotation() {
    constexpr double THETA_X = 1e-6;
    const SO3::AlgebraVector around_x_axis(THETA_X, 0., 0.);
    const linalg::Matrix3d rot_x = SO3::exp(around_x_axis).as_matrix();

    constexpr double TOL = 1e-10;
    util::assert_near(rot_x(1, 1), std::cos(THETA_X), TOL);
    util::assert_near(rot_x(1, 2), -std::sin(THETA_X), TOL);
    util::assert_near(rot_x(2, 1), std::sin(THETA_X), TOL);
    util::assert_near(rot_x(2, 2), std::cos(THETA_X), TOL);
    assert(rot_x(0, 0) == 1.0); assert(rot_x(0, 1) == 0.0);
    assert(rot_x(0, 2) == 0.0); assert(rot_x(1, 0) == 0.0);
    assert(rot_x(2, 0) == 0.0);
}

void exp_of_non_axis_aligned() {
    const SO3::AlgebraVector axis_angle(0.1, -0.75, 1.3);
    const linalg::Matrix3d skew = SO3::skew_matrix(axis_angle);
    linalg::Matrix3d nth_term = linalg::Matrix3d::identity();

    linalg::Matrix3d exp_numerical = linalg::Matrix3d::zero();
    constexpr int NUM_TERMS = 50;
    for (int i = 0; i < NUM_TERMS; ++i) {
        exp_numerical += nth_term;
        nth_term = (nth_term / (i + 1)) * skew;
    }

    const SO3 exp_analytic = SO3::exp(axis_angle);

    constexpr double TOL = 1e-10;
    assert_near(exp_analytic, SO3(exp_numerical), TOL);
}

void log_of_identity_is_zero() {
    const SO3::AlgebraVector should_be_zero = SO3().log();
    for (int i = 0; i < 3; ++i) {
        assert(should_be_zero(i) == 0.);
    }
}

void log_of_rotation_around_x() {
    constexpr double THETA_X = -M_PI / 7;
    const SO3::AlgebraVector around_x_axis(THETA_X, 0., 0.);
    const SO3 rot_x = make_rotation_about_x(THETA_X);
    const SO3::AlgebraVector should_be_orig = rot_x.log();

    constexpr double TOL = 1e-10;
    for (int i = 0; i < 3; ++i) {
        util::assert_near(around_x_axis(i), should_be_orig(i), TOL);
    }
}

void log_of_rotation_around_y() {
    constexpr double THETA_Y =  M_PI / 3;
    const SO3::AlgebraVector around_y_axis(0., THETA_Y, 0.);
    const SO3 rot_y = make_rotation_about_y(THETA_Y);
    const SO3::AlgebraVector should_be_orig = rot_y.log();

    constexpr double TOL = 1e-10;
    for (int i = 0; i < 3; ++i) {
        util::assert_near(around_y_axis(i), should_be_orig(i), TOL);
    }
}

void log_of_rotation_around_z() {
    constexpr double THETA_Z = - M_PI / 10;
    const SO3::AlgebraVector around_z_axis(0., 0., THETA_Z);
    const SO3 rot_z = make_rotation_about_z(THETA_Z);
    const SO3::AlgebraVector should_be_orig = rot_z.log();

    constexpr double TOL = 1e-10;
    for (int i = 0; i < 3; ++i) {
        util::assert_near(around_z_axis(i), should_be_orig(i), TOL);
    }
}

void log_of_small_rotation() {
    constexpr double THETA_Z = 1e-16;
    const SO3::AlgebraVector around_z_axis(0., 0., THETA_Z);
    const SO3 rot_z = make_rotation_about_z(THETA_Z);
    const SO3::AlgebraVector should_be_orig = rot_z.log();

    constexpr double TOL = 1e-16;
    for (int i = 0; i < 3; ++i) {
        util::assert_near(around_z_axis(i), should_be_orig(i), TOL);
    }
}

void log_of_rotation_by_pi() {
    const SO3::AlgebraVector axis_angle =
        SO3::AlgebraVector(1.2, 0.1, 0.7).normalized() * M_PI;
    const SO3 rot = SO3::exp(axis_angle);
    const SO3::AlgebraVector should_be_orig = rot.log();

    constexpr double TOL = 1e-10;
    for (int i = 0; i < 3; ++i) {
        util::assert_near(axis_angle(i), should_be_orig(i), TOL);
    }
}

void log_of_angle_somewhat_close_to_pi() {
    constexpr double THETA = M_PI - 1e-3;
    const SO3::AlgebraVector axis_angle =
        SO3::AlgebraVector(1.2, 0.1, 0.7).normalized() * THETA;
    const SO3 rot = SO3::exp(axis_angle);
    const SO3::AlgebraVector should_be_orig = rot.log();

    constexpr double TOL = 1e-10;
    for (int i = 0; i < 3; ++i) {
        util::assert_near(axis_angle(i), should_be_orig(i), TOL);
    }
}

void log_of_angle_very_close_to_pi() {
    constexpr double THETA = M_PI - 1e-16;
    const SO3::AlgebraVector axis_angle =
        SO3::AlgebraVector(1.2, 0.1, 0.7).normalized() * THETA;
    const SO3 rot = SO3::exp(axis_angle);
    const SO3::AlgebraVector should_be_orig = rot.log();

    constexpr double TOL = 1e-10;
    for (int i = 0; i < 3; ++i) {
        util::assert_near(axis_angle(i), should_be_orig(i), TOL);
    }
}

void differential_of_exp() {
    const SO3::AlgebraVector axis_angle(0.1, -0.75, 1.3);
    SO3::AlgebraTransformation diff_analytic;
    SO3::exp(axis_angle, &diff_analytic);
    constexpr double STEP_SIZE = 1e-6;
    const SO3::AlgebraTransformation diff_numerical =
        util::exp_diff_numerical<SO3>(axis_angle, STEP_SIZE);
    constexpr double TOL = 1e-6;
    for (int i  = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            util::assert_near(diff_analytic(i, j), diff_numerical(i, j), TOL);
        }
    }
}

void differential_of_log() {
    const SO3 rot = SO3::exp(SO3::AlgebraVector(0.25, -0.12, 0.33));
    SO3::AlgebraTransformation diff_analytic;
    rot.log(&diff_analytic);
    constexpr double STEP_SIZE = 1e-6;
    const SO3::AlgebraTransformation diff_numerical =
        util::log_diff_numerical<SO3>(rot, STEP_SIZE);
    constexpr double TOL = 1e-6;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            util::assert_near(diff_analytic(i, j), diff_numerical(i, j), TOL);
        }
    }
}

void rectified_identity_is_identity() {
    const SO3 should_be_identity = SO3().rectified();
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            if (i == j) {
                assert(should_be_identity.as_matrix()(i, j) == 1.);
            } else {
                assert(should_be_identity.as_matrix()(i, j) == 0.);
            }
        }
    }
}

void rectified_is_orthogonal() {
    constexpr double PERTURB = 1e-3;
    const SO3 almost_orthogonal =
        SO3(SO3::exp(SO3::AlgebraVector(0.01, 0.09, 0.99)).as_matrix() +
            linalg::Matrix3d::filled(PERTURB));
    const SO3 should_be_orthogonal = almost_orthogonal.rectified();

    constexpr double TOL = 1e-10;
    const linalg::Matrix3d should_be_identity =
        (should_be_orthogonal * should_be_orthogonal.inverse()).as_matrix();
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            if (i == j) {
                util::assert_near(should_be_identity(i, j), 1., TOL);
            } else {
                util::assert_near(should_be_identity(i, j), 0., TOL);
            }
        }
    }
    assert_near(almost_orthogonal, should_be_orthogonal, 3 * PERTURB);
}

}
}

int main() {
    std::cout << "exp_of_zero_is_identity..." << std::endl;
    sfm::liegroups::exp_of_zero_is_identity();

    std::cout << "exp_of_rotation_around_x..." << std::endl;
    sfm::liegroups::exp_of_rotation_around_x();

    std::cout << "exp_of_rotation_around_y..." << std::endl;
    sfm::liegroups::exp_of_rotation_around_y();

    std::cout << "exp_of_rotation_around_z..." << std::endl;
    sfm::liegroups::exp_of_rotation_around_z();

    std::cout << "exp_of_small_angle_rotation..." << std::endl;
    sfm::liegroups::exp_of_small_angle_rotation();

    std::cout << "exp_of_non_axis_aligned..." << std::endl;
    sfm::liegroups::exp_of_non_axis_aligned();

    std::cout << "log_of_identity_is_zero..." << std::endl;
    sfm::liegroups::log_of_identity_is_zero();

    std::cout << "log_of_rotation_around_x..." << std::endl;
    sfm::liegroups::log_of_rotation_around_x();

    std::cout << "log_of_rotation_around_y..." << std::endl;
    sfm::liegroups::log_of_rotation_around_y();

    std::cout << "log_of_rotation_around_z..." << std::endl;
    sfm::liegroups::log_of_rotation_around_z();

    std::cout << "log_of_small_rotation..." << std::endl;
    sfm::liegroups::log_of_small_rotation();

    std::cout << "log_of_rotation_by_pi..." << std::endl;
    sfm::liegroups::log_of_rotation_by_pi();

    std::cout << "log_of_angle_somewhat_close_to_pi..." << std::endl;
    sfm::liegroups::log_of_angle_somewhat_close_to_pi();

    std::cout << "log_of_angle_very_close_to_pi" << std::endl;
    sfm::liegroups::log_of_angle_very_close_to_pi();

    std::cout << "differential_of_exp..." << std::endl;
    sfm::liegroups::differential_of_exp();

    std::cout << "differential_of_log..." << std::endl;
    sfm::liegroups::differential_of_log();

    std::cout << "rectified_identity_is_identity..." << std::endl;
    sfm::liegroups::rectified_identity_is_identity();

    std::cout << "rectified_is_orthogonal..." << std::endl;
    sfm::liegroups::rectified_is_orthogonal();

    std::cout << "All tests passed!" << std::endl;
    return 0;
}
