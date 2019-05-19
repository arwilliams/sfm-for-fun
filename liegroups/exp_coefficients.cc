#include "liegroups/exp_coefficients.hh"

#include <cmath>

namespace sfm {
namespace liegroups {
namespace {

constexpr double SMALL_THETA = 1e-8;

double compute_a(const double theta) {
	if (theta < SMALL_THETA) {
		const double theta_sq = theta * theta;
		return 1. - theta_sq / 6. + theta_sq * theta_sq / 120.;
	}
	return std::sin(theta) / theta;
}

double compute_b(const double theta) {
	if (theta < SMALL_THETA) {
		const double theta_sq = theta * theta;
		return 0.5 - theta_sq / 24. + theta_sq * theta_sq / 720.;
	}
	return (1 - std::cos(theta)) / (theta * theta);
}

double compute_c(const double theta) {
    if (theta < SMALL_THETA) {
        const double theta_sq = theta * theta;
        return 1. / 6. - theta_sq / 120. + theta_sq * theta_sq / 5040.;
    }
    return (1. - compute_a(theta)) / (theta * theta);
}

}

ExpCoefficients compute_exp_coefficients(const double theta) {
    ExpCoefficients coeffs;
    coeffs.a = compute_a(theta);
    coeffs.b = compute_b(theta);
    coeffs.c = compute_c(theta);
    
    const double cos_theta = std::cos(theta);
    if (cos_theta > 0.999856) {
        coeffs.e = (coeffs.b - 2 * coeffs.c) / (2. * coeffs.a);
    } else {
        coeffs.e = (coeffs.b - 0.5 * coeffs.a) / (1. - cos_theta);
    }

    return coeffs;
}

}
}

