#pragma once

#include "linalg/matrix.hh"

#include <cassert>
#include <cmath>

namespace sfm {
namespace linalg {

template <class Derived>
bool cholesky_in_place(MatrixBase<Derived> &A) {
    constexpr double SMALL_VALUE = 1e-10;
    ASSERT_SQUARE(A);
    for (int i = 0; i < A.rows(); ++i) {
        for (int j = 0; j < A.cols(); ++j) {
            if (j > i) {
                A(i, j) = 0;
                continue;
            }
            double sum = 0.;
            for (int k = 0; k < j; ++k) {
                if (i == j) {
                    sum += A(j, k) * A(j, k);
                } else {
                    sum += A(i, k) * A(j, k);
                }
            }
            A(i, j) -= sum;
            if (i == j) {
				if (A(i, j) < SMALL_VALUE) {
                    return false;
                }
                A(i, j) = std::sqrt(A(i, j));
            } else {
				if (A(j, j) < SMALL_VALUE) {
                    return false;
                }
                A(i, j) /= A(j, j);
            }
        }
    }
    return true;
}

// solve Ly = b for y, given lower-triangular L
template <class Derived, class OtherDerived>
auto solve_forward_sub(
        const MatrixBase<Derived> &L,
        const MatrixBase<OtherDerived> &b) {
    ASSERT_SQUARE(L);
    static constexpr int DIM = traits<Derived>::rows;
    ASSERT_IS_VECTOR<DIM>(b);
    ASSERT_SAME_SCALAR_TYPE(L, b);
    constexpr double SMALL_VALUE = 1e-10;

    Vector<typename traits<Derived>::scalar_t, DIM> y;
	for (int i = 0; i < y.rows(); ++i) {
		y(i) = b(i, 0);
		for (int k = 0; k < i; ++k) {
			y(i) -= L(i, k) * y(k);
		}
		assert(L(i, i) > SMALL_VALUE);
		y(i) /= L(i, i);
	}
	return y;
}

// solve Ux = y, given upper-triangular U
template <class Derived, class OtherDerived>
auto solve_back_sub(const MatrixBase<Derived> &U,
                    const MatrixBase<OtherDerived> &y) {
    ASSERT_SQUARE(U);
    static constexpr int DIM = traits<Derived>::rows;
    ASSERT_IS_VECTOR<DIM>(y);
    ASSERT_SAME_SCALAR_TYPE(U, y);
    auto x = Vector<typename traits<Derived>::scalar_t, DIM>::zero();
    constexpr double SMALL_VALUE = 1e-10;

	for (int i = y.rows() - 1; i >= 0; --i) {
		x(i) = y(i);
		for (int k = i + 1; k < x.rows(); ++k) {
			x(i) -= U(i, k) * x(k);
		}
		assert(U(i, i) > SMALL_VALUE);
		x(i) /= U(i, i);
	}
	return x;
}

enum class CholeskyStatus : uint8_t {
    DECOMPOSITION_FAILED,
    OKAY,
};

template <class Derived>
class CholeskyDecomposition {
 public:
    static constexpr int DIM = traits<Derived>::rows;

 	explicit CholeskyDecomposition(const MatrixBase<Derived> &A) {
        ASSERT_SQUARE(A);
        L_ = A;
        if (cholesky_in_place(L_)) {
            status_ = CholeskyStatus::OKAY;
        } else {
            status_ = CholeskyStatus::DECOMPOSITION_FAILED;
        }
    }

    CholeskyStatus status() const { return status_; }

    template <class OtherDerived>
	auto solve(const MatrixBase<OtherDerived> &b) const {
        assert(status_ == CholeskyStatus::OKAY);
		// Ly = b
		const auto y = solve_forward_sub(L_, b);
		// L^Tx = y
		const auto x = solve_back_sub(L_.transpose(), y);
		// => Ax = b
		return x;
	}

    const auto &L() const { return L_; }
 private:
    Matrix<typename traits<Derived>::scalar_t, DIM, DIM> L_;
    CholeskyStatus status_;
};

}  // namespace linalg
}  // namespace sfm
