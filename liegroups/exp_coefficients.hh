#pragma once

namespace sfm {
namespace liegroups {

//
// Inspired by the notation used at ethaneade.com.
//
struct ExpCoefficients {
    // sin(t) / t
    double a;

    // (1 - cos(t)) / t^2
    double b;

    // (1 - a(t)) / t^2
    double c;

    // (b(t) - 2c(t)) / (2a(t))
    // = (b(t) - 0.5a(t)) / (1 - cos(t))
    double e;
};

ExpCoefficients compute_exp_coefficients(const double theta);

}
}
