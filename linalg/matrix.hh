#pragma once

#include "linalg/forward_declarations.hh"

#include <array>
#include <cassert>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <type_traits>

namespace sfm {
namespace linalg {

template <class Derived>
struct traits<MatrixBase<Derived>> {
	typedef traits<Derived> traits_t;
	typedef typename traits_t::scalar_t scalar_t;
	typedef typename traits_t::transpose_t transpose_t;
	static constexpr int rows = traits_t::rows;
	static constexpr int cols = traits_t::cols;
};

template <typename Scalar, int Rows, int Cols>
struct traits<Matrix<Scalar, Rows, Cols>> {
	typedef Matrix<Scalar, Rows, Cols> matrix_t;
	typedef TransposeMatrix<matrix_t> transpose_t;
	typedef Scalar scalar_t;
	static constexpr int rows = Rows;
	static constexpr int cols = Cols;
};

template <class Derived>
struct traits<TransposeMatrix<Derived>> {
	typedef typename traits<Derived>::scalar_t scalar_t;
	typedef const Derived &transpose_t;

	static constexpr int rows = traits<Derived>::cols;
	static constexpr int cols = traits<Derived>::rows;
};

template <class Derived, int Rows, int Cols>
struct traits<Block<Derived, Rows, Cols>> {
	typedef typename traits<Derived>::scalar_t scalar_t;
	typedef Block<Derived, Rows, Cols> matrix_t;
	static constexpr int rows = Rows;
	static constexpr int cols = Cols;
    typedef TransposeMatrix<matrix_t> transpose_t;
};

template <class Derived, class OtherDerived, class BinaryOp>
struct traits<CmptWiseBinaryMatrixOpResult<Derived, OtherDerived, BinaryOp>> {
	typedef typename traits<Derived>::scalar_t scalar_t;
	typedef CmptWiseBinaryMatrixOpResult<Derived, OtherDerived, BinaryOp> matrix_t;

	static constexpr int rows = traits<Derived>::rows > traits<OtherDerived>::rows ?
									   	traits<Derived>::rows
									   : traits<OtherDerived>::rows;
	static constexpr int cols = traits<Derived>::cols > traits<OtherDerived>::cols ?
									   	traits<Derived>::cols
									   : traits<OtherDerived>::cols;

	typedef TransposeMatrix<matrix_t> transpose_t;
};

template <class Derived, class Scalar, class BinaryOp>
struct traits<CmptWiseMatrixScalarOpResult<Derived, Scalar, BinaryOp>> {
	typedef typename traits<Derived>::scalar_t scalar_t;
	typedef CmptWiseMatrixScalarOpResult<Derived, Scalar, BinaryOp> matrix_t;

	static constexpr int rows = traits<Derived>::rows;
	static constexpr int cols = traits<Derived>::cols;

	typedef TransposeMatrix<matrix_t> transpose_t;
};

template <class Derived, class UnaryOp>
struct traits<CmptWiseUnaryMatrixOpResult<Derived, UnaryOp>> {
	typedef typename traits<Derived>::scalar_t scalar_t;
	typedef CmptWiseUnaryMatrixOpResult<Derived, UnaryOp> matrix_t;

	static constexpr int rows = traits<Derived>::rows;
	static constexpr int cols = traits<Derived>::cols;

	typedef TransposeMatrix<matrix_t> transpose_t;
};

template <class Derived, class OtherDerived>
struct traits<MatrixProductResult<Derived, OtherDerived>> {
	typedef typename traits<Derived>::scalar_t scalar_t;
	typedef MatrixProductResult<Derived, OtherDerived> matrix_t;

	static constexpr int rows = traits<Derived>::rows;
	static constexpr int cols = traits<OtherDerived>::cols;

	typedef TransposeMatrix<matrix_t> transpose_t;
};

template <class Derived>
void assert_in_bounds(int i, int j, const MatrixBase<Derived> &a) {
    if (!(i >= 0 && i < a.rows() && j >= 0 && j < a.cols())) {
        std::stringstream ss;
        ss << "Invalid index: (" << i << ", " << j << ") into matrix with dims ("
           << a.rows() << ", " << a.cols() << ")";
        throw std::runtime_error(ss.str());
    }
}

template <class Derived>
void assert_in_bounds(int i, const MatrixBase<Derived> &a) {
    static_assert(traits<Derived>::rows == 1 || traits<Derived>::cols == 1);
    if (!(i >= 0 && i < a.rows() * a.cols())) {
        std::stringstream ss;
        ss << "Invalid index: " << i << " into vector with dims ("
           << a.rows() << ", " << a.cols() << ")";
        throw std::runtime_error(ss.str());
    }
}

template <class Derived, class OtherDerived>
void ASSERT_SAME_SHAPE(const MatrixBase<Derived> &a, const MatrixBase<OtherDerived> &b) {
	static_assert(traits<Derived>::rows == traits<OtherDerived>::rows,
				  "Matrix dimension mismatch (rows)");
	static_assert(traits<Derived>::cols == traits<OtherDerived>::cols,
				  "Matrix dimension mismatch (cols)");
}

template <class Derived, class OtherDerived>
void ASSERT_SAME_SCALAR_TYPE(const MatrixBase<Derived> &a, const MatrixBase<OtherDerived> &b) {
	static_assert(std::is_same<typename traits<Derived>::scalar_t, typename traits<OtherDerived>::scalar_t>::value,
			      "Matrix scalar type mismatch");
}

template <class Derived, class OtherDerived>
void ASSERT_MULTIPLIABLE(const MatrixBase<Derived> &a, const MatrixBase<OtherDerived> &b) {
	static_assert(traits<Derived>::cols == traits<OtherDerived>::rows,
				  "Matrix mutliplication dimension mismatch");
}

template <class Derived>
void ASSERT_SQUARE(const MatrixBase<Derived> &a) {
	static_assert(traits<Derived>::rows == traits<Derived>::cols,
				  "Non-square matrix");
}

template <int Cols, class Derived>
void ASSERT_NUM_COLS(const MatrixBase<Derived> &a) {
    static_assert(traits<Derived>::cols == Cols, "Wrong num cols");
}


template <int Rows, class Derived>
void ASSERT_NUM_ROWS(const MatrixBase<Derived> &a) {
    static_assert(traits<Derived>::rows == Rows, "Wrong num rows");
}

template <int Dim, class Derived>
void ASSERT_IS_VECTOR(const MatrixBase<Derived> &a) {
    ASSERT_NUM_COLS<1>(a);
    ASSERT_NUM_ROWS<Dim>(a);
}

template <class Derived>
class MatrixBase {
 public:
 	Derived &derived() { return static_cast<Derived &>(*this); }
	const Derived &derived() const { return static_cast<const Derived &>(*this); }

	typedef typename traits<Derived>::scalar_t scalar_t;
	typedef typename traits<Derived>::transpose_t transpose_t;

	static constexpr int rows() { return traits<Derived>::rows; }
	static constexpr int cols() { return traits<Derived>::cols; }

	scalar_t operator()(int i, int j) const { return derived()(i, j); }
	scalar_t &operator()(int i, int j) { return derived()(i, j); }

    scalar_t operator()(int i) const {
        static_assert(traits<Derived>::rows == 1 || traits<Derived>::cols == 1);
        return derived()(i);
    }

    scalar_t &operator()(int i) {
        static_assert(traits<Derived>::rows == 1 || traits<Derived>::cols == 1);
        return derived()(i);
    }

	transpose_t transpose() const { return derived().transpose(); }

    template<int Rows, int Cols>
    auto block(int i0, int j0) const { return derived().block<Rows, Cols>(i0, j0); }

    auto row(int i) const { return derived().row(i); }
    auto row(int i) { return derived().row(i); }

    auto col(int j) const { return derived().col(j); }
    auto col(int j) { return derived().col(j); }

    MatrixBase &operator=(const MatrixBase<Derived> &other) {
        this->safe_assign(other);
        return *this;
    }

	template <class OtherDerived>
	MatrixBase &operator=(const MatrixBase<OtherDerived> &other) {
		this->safe_assign(other);
		return *this;
	}

	template <class OtherDerived>
	MatrixBase &operator+=(const MatrixBase<OtherDerived> &other) {
        this->safe_assign(*this + other);
        return *this;
    }

	template <class OtherDerived>
	MatrixBase &operator-=(const MatrixBase<OtherDerived> &other) {
        this->safe_assign(*this - other);
        return *this;
    }

	template <typename Scalar>
	MatrixBase &operator*=(Scalar s) {
        this->safe_assign(*this * s);
        return *this;
    }

	template <typename Scalar>
	MatrixBase &operator/=(Scalar s) {
        this->safe_assign(*this / s);
        return *this;
    }

    template <typename Scalar>
    void fill(Scalar s) {
        for (int i = 0; i < rows(); ++i) {
            for (int j = 0; j < cols(); ++j) {
                (*this)(i, j) = s;
            }
        }
    }

    static Derived zero() {
        Derived a;
        a.fill(scalar_t(0));
        return a;
    }

    static Derived identity() {
        Derived a;
        ASSERT_SQUARE(a);
        for (int i = 0; i < a.rows(); ++i) {
            for (int j = 0; j < a.cols(); ++j) {
                a(i, j) = static_cast<int>(i == j);
            }
        }
        return a;
    }

    static Derived filled(const scalar_t s) {
        Derived a;
        a.fill(s);
        return a;
    }

    scalar_t squared_norm() const {
        scalar_t nsq(0);
        for (int i = 0; i < rows(); ++i) {
            for (int j = 0; j < cols(); ++j) {
                nsq += (*this)(i, j) * (*this)(i, j);
            }
        }
        return nsq;
    }

    scalar_t norm() const {
        return std::sqrt(squared_norm());
    }

    scalar_t trace() const {
        ASSERT_SQUARE(*this);
        scalar_t tr = 0;
        for (int i = 0; i < rows(); ++i) {
            tr += (*this)(i, i);
        }
        return tr;
    }

    void normalize() {
        *this /= norm();
    }

    Derived normalized() const {
        Derived out = this->derived();
        out.normalize();
        return out;
    }

    template <class OtherDerived>
    scalar_t dot(const OtherDerived &other) const {
        ASSERT_IS_VECTOR<rows()>(*this);
        ASSERT_IS_VECTOR<rows()>(other);
        scalar_t result = 0;
        for (int i = 0; i < rows(); ++i) {
            result += (*this)(i) * other(i);
        }
        return result;
    }

 protected:
 	template <class OtherDerived>
	void assign(const MatrixBase<OtherDerived> &other) {
		for (int i = 0; i < rows(); ++i) {
			for (int j = 0; j < cols(); ++j) {
				(*this)(i, j) = other(i, j);
			}
		}
	}

 	template <class OtherDerived>
	void safe_assign(const MatrixBase<OtherDerived> &other) {
		ASSERT_SAME_SHAPE(*this, other);
		this->assign(other);
	}
};

template <typename Scalar, int Rows, int Cols>
class Matrix : public MatrixBase<Matrix<Scalar, Rows, Cols>> {
 public:
 	typedef Matrix<Scalar, Rows, Cols> self_t;
 	typedef typename traits<self_t>::scalar_t scalar_t;
 	typedef typename traits<self_t>::transpose_t transpose_t;

	Matrix() = default;

	scalar_t operator()(int i, int j) const {
        assert_in_bounds(i, j, *this);
		return data_[Cols * i + j];
	}

	scalar_t &operator()(int i, int j) {
        assert_in_bounds(i, j, *this);
		return data_[Cols * i + j];
	}

    scalar_t &operator()(int i) {
        assert_in_bounds(i, *this);
        return data_[i];
    }

    scalar_t operator()(int i) const {
        assert_in_bounds(i, *this);
        return data_[i];
    }

	auto transpose() const { return transpose_t(*this); }
    
    template <int BlockRows, int BlockCols>
    auto block(int i0, int j0) const {
        return Block<self_t, BlockRows, BlockCols>(*this, i0, j0);
    }

    template <int BlockRows, int BlockCols>
    auto block(int i0, int j0) {
        return Block<self_t, BlockRows, BlockCols>(*this, i0, j0);
    }

    auto col(int j) const { return Block<self_t, Rows, 1>(*this, 0, j); }
    auto col(int j) { return Block<self_t, Rows, 1>(*this, 0, j); }

    auto row(int i) const { return Block<self_t, 1, Cols>(*this, i, 0); }
    auto row(int i) { return Block<self_t, 1, Cols>(*this, i, 0); }

	template <class OtherDerived>
	Matrix(const MatrixBase<OtherDerived> &other) {
		this->safe_assign(other);
	}

    // Specialized constructors for vector types.
    explicit Matrix(const scalar_t x) : data_{x} {
        ASSERT_IS_VECTOR<1>(*this);
    }

    Matrix(const scalar_t x, const scalar_t y)
        : data_{x, y} {
        ASSERT_IS_VECTOR<2>(*this);
    }

    Matrix(const scalar_t x, const scalar_t y, const scalar_t z)
        : data_{x, y, z} {
        ASSERT_IS_VECTOR<3>(*this);
    }

    static Matrix<Scalar, Rows, Cols> unit_vec(int i) {
        Matrix<Scalar, Rows, Cols> x =
            Matrix<Scalar, Rows, Cols>::zero();
        ASSERT_IS_VECTOR<Rows>(x);
        x(i) = 1.0;
        return x;
    }

 private:
 	std::array<Scalar, Rows * Cols> data_;
};

template <typename Scalar, int Rows>
using Vector = Matrix<Scalar, Rows, 1>;

template <class Derived>
class TransposeMatrix : public MatrixBase<TransposeMatrix<Derived>> {
 public:
 	TransposeMatrix(const Derived &mat) : mat_(mat) {}

	typedef TransposeMatrix<Derived> self_t;
	typedef typename traits<self_t>::scalar_t scalar_t;
	typedef typename traits<self_t>::transpose_t transpose_t;

	scalar_t operator()(int i, int j) const { return mat_(j, i); }

	transpose_t transpose() const { return mat_; }
 private:
 	const Derived &mat_;
};

template <class Derived, class OtherDerived, class BinaryOp>
class CmptWiseBinaryMatrixOpResult : public MatrixBase<CmptWiseBinaryMatrixOpResult<Derived, OtherDerived, BinaryOp>> {
 public:
 	CmptWiseBinaryMatrixOpResult(const Derived &a, const OtherDerived &b) : a_(a), b_(b) {
		ASSERT_SAME_SHAPE(a, b);
		ASSERT_SAME_SCALAR_TYPE(a, b);
	}

	const BinaryOp op = BinaryOp();

	typedef CmptWiseBinaryMatrixOpResult<Derived, OtherDerived, BinaryOp> self_t;
	typedef typename traits<self_t>::scalar_t scalar_t;
	typedef typename traits<self_t>::transpose_t transpose_t;

	scalar_t operator()(int i, int j) const { return op(a_(i, j),  b_(i, j)); }
 private:
 	const Derived &a_;
	const OtherDerived &b_;
};

template <class Derived, typename Scalar, class BinaryOp>
class CmptWiseMatrixScalarOpResult : public MatrixBase<CmptWiseMatrixScalarOpResult<Derived, Scalar, BinaryOp>> {
 public:
 	CmptWiseMatrixScalarOpResult(const Derived &a, Scalar s) : a_(a), s_(s) {}

	const BinaryOp op = BinaryOp();

	typedef CmptWiseMatrixScalarOpResult<Derived, Scalar, BinaryOp> self_t;
	typedef typename traits<self_t>::scalar_t scalar_t;

	scalar_t operator()(int i, int j) const { return op(a_(i, j), static_cast<scalar_t>(s_)); }
 private:
	const Derived &a_;
 	const Scalar s_;
};

template <class Derived, class UnaryOp>
class CmptWiseUnaryMatrixOpResult : public MatrixBase<CmptWiseUnaryMatrixOpResult<Derived, UnaryOp>> {
 public:
 	CmptWiseUnaryMatrixOpResult(const Derived &a) : a_(a) {}

	const UnaryOp op = UnaryOp();

	typedef CmptWiseUnaryMatrixOpResult<Derived, UnaryOp> self_t;
	typedef typename traits<self_t>::scalar_t scalar_t;
	typedef typename traits<self_t>::transpose_t transpose_t;

	scalar_t operator()(int i, int j) const { return op(a_(i, j)); }
 private:
 	const Derived &a_;
};

template <class Derived, class OtherDerived>
class MatrixProductResult : public MatrixBase<MatrixProductResult<Derived, OtherDerived>> {
 public:
 	MatrixProductResult(const Derived &a, const OtherDerived &b) : a_(a), b_(b) {
		ASSERT_SAME_SCALAR_TYPE(a, b);
		ASSERT_MULTIPLIABLE(a, b);
	}

	typedef MatrixProductResult<Derived, OtherDerived> self_t;
	typedef typename traits<self_t>::scalar_t scalar_t;

	scalar_t operator()(int i, int j) const {
		scalar_t sum = 0;
		for (int k = 0; k < a_.cols(); ++k) {
			sum += a_(i, k) * b_(k, j);
		}
		return sum;
	}
 private:
	const Derived &a_;
 	const OtherDerived &b_;
};

template <class Derived, int Rows, int Cols>
class Block : public MatrixBase<Block<Derived, Rows, Cols>> {
 public:
    Block(const Derived &a, int i0, int j0)
        : const_a_(a), i0_(i0), j0_(j0) {}
    Block(Derived &a, int i0, int j0)
        : const_a_(a), a_(&a), i0_(i0), j0_(j0) {}

     typedef Block<Derived, Rows, Cols> self_t;
     typedef typename traits<self_t>::scalar_t scalar_t;

     using MatrixBase<Block<Derived, Rows, Cols>>::operator=;

     Block &operator=(const Block<Derived, Rows, Cols> &other) {
         return static_cast<Block &>(MatrixBase<Block<Derived, Rows, Cols>>::operator=(other));
     }

     scalar_t operator()(int i, int j) const {
         return const_a_(i0_ + i, j0_ + j);
     }

     scalar_t &operator()(int i, int j) {
         return (*a_)(i0_ + i, j0_ + j);
     }

     const Derived &const_a_;
     Derived *a_;
     int i0_;
     int j0_;
};

template <class Derived, class OtherDerived>
using MatrixAddResult = CmptWiseBinaryMatrixOpResult<Derived, OtherDerived, std::plus<typename traits<Derived>::scalar_t>>;

template <class Derived, class OtherDerived>
MatrixAddResult<Derived, OtherDerived> operator+(
    const MatrixBase<Derived> &a, const MatrixBase<OtherDerived> &b) {
	return MatrixAddResult<Derived, OtherDerived>(a.derived(), b.derived());
}

template <class Derived, class OtherDerived>
using MatrixSubtractResult = CmptWiseBinaryMatrixOpResult<Derived, OtherDerived, std::minus<typename traits<Derived>::scalar_t>>;

template <class Derived, class OtherDerived>
MatrixSubtractResult<Derived, OtherDerived> operator-(const MatrixBase<Derived> &a, const MatrixBase<OtherDerived> &b) {
	return MatrixSubtractResult<Derived, OtherDerived>(a.derived(), b.derived());
}

template <class Derived, typename Scalar>
using MatrixScalarMultResult = CmptWiseMatrixScalarOpResult<Derived, Scalar, std::multiplies<typename traits<Derived>::scalar_t>>;

template <class Derived, typename Scalar,
          typename = std::enable_if_t<std::is_arithmetic<Scalar>::value>>
MatrixScalarMultResult<Derived, Scalar> operator*(const MatrixBase<Derived> &a, const Scalar s) {
	return MatrixScalarMultResult<Derived, Scalar>(a.derived(), s);
}

template <class Derived, typename Scalar,
          typename = std::enable_if_t<std::is_arithmetic<Scalar>::value>>
MatrixScalarMultResult<Derived, Scalar> operator*(const Scalar s, const MatrixBase<Derived> &a) {
	return MatrixScalarMultResult<Derived, Scalar>(a.derived(), s);
}

template <class Derived, typename Scalar>
using MatrixScalarDivResult = CmptWiseMatrixScalarOpResult<Derived, Scalar, std::divides<typename traits<Derived>::scalar_t>>;

template <class Derived, typename Scalar>
MatrixScalarDivResult<Derived, Scalar> operator/(const MatrixBase<Derived> &a, const Scalar s) {
	return MatrixScalarDivResult<Derived, Scalar>(a.derived(), s);
}

template <class Derived, class OtherDerived>
MatrixProductResult<Derived, OtherDerived> operator*(
    const MatrixBase<Derived> &a, const MatrixBase<OtherDerived> &b) {
    return MatrixProductResult<Derived, OtherDerived>(a.derived(), b.derived());
}

template <class Derived>
using MatrixNegateResult = CmptWiseUnaryMatrixOpResult<Derived, std::negate<typename traits<Derived>::scalar_t>>;

template <class Derived>
MatrixNegateResult<Derived> operator-(const MatrixBase<Derived> &a) {
    return MatrixNegateResult<Derived>(a.derived());
}

using Matrix3d = Matrix<double, 3, 3>;
using Vector3d = Matrix<double, 3, 1>;

using Matrix4d = Matrix<double, 4, 4>;
using Vector4d = Matrix<double, 4, 1>;

}  // namespace linalg
}  // namespace sfm
