#include "linalg/cholesky.hh"

#include <cassert>
#include <iostream>

namespace sfm {
namespace linalg {

void cholesky_one_by_one_is_square_root() {
    constexpr double VAL = 2.0;
    Matrix<double, 1, 1> x(VAL);
    const CholeskyDecomposition<decltype(x)> llt(x);
    assert(llt.status() == CholeskyStatus::OKAY);
    assert(llt.L()(0, 0) == std::sqrt(VAL));
}

void cholesky_succeeds_on_valid_input_int_matrix() {
    Matrix<int, 3, 3> m;
    m(0, 0) = 4;   m(0, 1) = 12;  m(0, 2) = -16;
    m(1, 0) = 12;  m(1, 1) = 37;  m(1, 2) = -43;
    m(2, 0) = -16; m(2, 1) = -43; m(2, 2) = 98;

    const CholeskyDecomposition<decltype(m)> llt(m);
    assert(llt.status() == CholeskyStatus::OKAY);
    const auto &L = llt.L();
    assert(L(0, 0) == 2);  assert(L(0, 1) == 0); assert(L(0, 2) == 0);
    assert(L(1, 0) == 6);  assert(L(1, 1) == 1); assert(L(1, 2) == 0);
    assert(L(2, 0) == -8); assert(L(2, 1) == 5); assert(L(2, 2) == 3);
}

void cholesky_succeeds_on_valid_input_double_matrix() {
    Matrix<double, 3, 3> m;
    m(0, 0) = 2.0;  m(0, 1) = -1.0; m(0, 2) = 0.0;
    m(1, 0) = -1.0; m(1, 1) = 2.0;  m(1, 2) = -1.0;
    m(2, 0) = 0.0;  m(2, 1) = -1.0; m(2, 0) = 2.0;

    const CholeskyDecomposition<decltype(m)> llt(m);
}

void cholesky_fails_on_negative_eigenvalue() {
    Matrix<int, 3, 3> m;
    m(0, 0) = 4;   m(0, 1) = 12;  m(0, 2) = -16;
    m(1, 0) = 12;  m(1, 1) = 37;  m(1, 2) = -43;
    m(2, 0) = -16; m(2, 1) = -43; m(2, 2) = -98;

    const CholeskyDecomposition<decltype(m)> llt(m);
    assert(llt.status() == CholeskyStatus::DECOMPOSITION_FAILED);
}

void cholesky_fails_on_positive_semidefinite() {
    const Matrix<double, 2, 2> m = Matrix<double, 2, 2>::filled(1.0);

    const CholeskyDecomposition<Matrix<double, 2, 2>> llt(m);
    assert(llt.status() == CholeskyStatus::DECOMPOSITION_FAILED);
}

}
}

int main() {
    std::cout << "cholesky_one_by_one_is_square_root..." << std::endl;
    sfm::linalg::cholesky_one_by_one_is_square_root();

    std::cout << "cholesky_succeeds_on_valid_input_int_matrix..." << std::endl;
    sfm::linalg::cholesky_succeeds_on_valid_input_int_matrix();

    std::cout << "cholesky_succeeds_on_valid_input_double_matrix..." << std::endl;
    sfm::linalg::cholesky_succeeds_on_valid_input_double_matrix();

    std::cout << "cholesky_fails_on_negative_eigenvalue..." << std::endl;
    sfm::linalg::cholesky_fails_on_negative_eigenvalue();

    std::cout << "cholesky_fails_on_positive_semidefinite..." << std::endl;
    sfm::linalg::cholesky_fails_on_positive_semidefinite();

    std::cout << "All tests passed!" << std::endl;
    return 0;
}
