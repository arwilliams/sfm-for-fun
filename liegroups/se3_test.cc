#include "liegroups/se3.hh"
#include "linalg/io.hh"

#include "util/test_helpers.hh"

#include <cassert>
#include <iostream>

namespace sfm {
namespace liegroups {
namespace {

void assert_near(const SE3 &x, const SE3 &y, const double tolerance) {
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            util::assert_near(x.as_matrix()(i, j), y.as_matrix()(i, j), tolerance);
        }
    }
}

void exp_of_zero_is_identity() {
    const SE3 should_be_identity = SE3::exp(SE3::AlgebraVector::zero());
    constexpr double TOL = 1e-20;
    assert_near(should_be_identity, SE3(), TOL);
}

void log_of_identity_is_zero() {
    const SE3::AlgebraVector should_be_zero = SE3().log();
    for (int i = 0; i < SE3::DOF; ++i) {
        assert(should_be_zero(i) == 0.0);
    }
}

void exp_matches_numerical() {
    const SO3::AlgebraVector w(0.12, -0.01, 0.79);
    const linalg::Vector3d u(12.3, -99.1, 56.3);
    SE3::AlgebraVector alg_vec;
    SE3::rotational_component(alg_vec) = w;
    SE3::translational_component(alg_vec) = u;

    linalg::Matrix4d alg_mat = linalg::Matrix4d::zero();
    alg_mat.block<3, 3>(0, 0) = SO3::adjoint(w);
    alg_mat.block<3, 1>(0, 3) = u;

    linalg::Matrix4d nth_term = linalg::Matrix4d::identity();

    linalg::Matrix4d exp_numerical_mat = linalg::Matrix4d::zero();
    constexpr int NUM_TERMS = 50;
    for (int i = 0; i < NUM_TERMS; ++i) {
        exp_numerical_mat += nth_term;
        nth_term = (nth_term / (i + 1)) * alg_mat;
    }

    const SE3 exp_numerical = SE3(SO3(exp_numerical_mat.block<3, 3>(0, 0)),
                                  exp_numerical_mat.block<3, 1>(0, 3));
                                
    const SE3 exp_analytic = SE3::exp(alg_vec);

    constexpr double TOL = 1e-10;
    assert_near(exp_analytic, exp_numerical, TOL);
}

void log_of_exp_is_identity() {
    const SO3::AlgebraVector w(0.12, -0.01, 0.79);
    const linalg::Vector3d u(12.3, -99.1, 56.3);
    SE3::AlgebraVector wu;
    SE3::rotational_component(wu) = w;
    SE3::translational_component(wu) = u;

    const SE3::AlgebraVector should_be_wu = SE3::exp(wu).log();
    constexpr double TOL = 1e-10;
    for (int i = 0; i < SE3::DOF; ++i) {
        util::assert_near(should_be_wu(i), wu(i), TOL);
    }
}

void differential_of_exp() {
    const SO3::AlgebraVector w(0.1, -0.75, 1.3);
    const linalg::Vector3d u(9.9, -100.0, 202.7);
    SE3::AlgebraVector wu;
    SE3::rotational_component(wu) = w;
    SE3::translational_component(wu) = u;

    SE3::AlgebraTransformation diff_analytic;
    SE3::exp(wu, &diff_analytic);
    constexpr double STEP_SIZE = 1e-8;
    const SE3::AlgebraTransformation diff_numerical =
        util::exp_diff_numerical<SE3>(wu, STEP_SIZE);
    constexpr double TOL = 1e-4;
    for (int i  = 0; i < SE3::DOF; ++i) {
        for (int j = 0; j < SE3::DOF; ++j) {
            util::assert_near(diff_analytic(i, j), diff_numerical(i, j), TOL);
        }
    }
}

void differential_of_log() {
    const SE3 x = SE3(SO3::exp(SO3::AlgebraVector(0.1, -0.75, 1.3)),
                      linalg::Vector3d(9.9, -100.0, 202.7));

    SE3::AlgebraTransformation diff_analytic;
    x.log(&diff_analytic);
    constexpr double STEP_SIZE = 1e-6;
    const SE3::AlgebraTransformation diff_numerical =
        util::log_diff_numerical<SE3>(x, STEP_SIZE);
    constexpr double TOL = 1e-4;
    for (int i  = 0; i < SE3::DOF; ++i) {
        for (int j = 0; j < SE3::DOF; ++j) {
            util::assert_near(diff_analytic(i, j), diff_numerical(i, j), TOL);
        }
    }
}

}
}
}

int main() {
    std::cout << "exp_of_zero_is_identity..." << std::endl;
    sfm::liegroups::exp_of_zero_is_identity();

    std::cout << "log_of_identity_is_zero..." << std::endl;
    sfm::liegroups::log_of_identity_is_zero();

    std::cout << "exp_matches_numerical..." << std::endl;
    sfm::liegroups::exp_matches_numerical();

    std::cout << "log_of_exp_is_identity..." << std::endl;
    sfm::liegroups::log_of_exp_is_identity();

    std::cout << "differential_of_exp..." << std::endl;
    sfm::liegroups::differential_of_exp();

    std::cout << "differential_of_log..." << std::endl;
    sfm::liegroups::differential_of_log();

    std::cout << "All tests passed!" << std::endl;
    return 0;
}
