#pragma once

#include "forward_declarations.hh"

#include <cassert>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <type_traits>

namespace sfm {

template <class Derived>
struct traits<MatrixBase<Derived>> {
	typedef traits<Derived> traits_t;
	typedef typename traits_t::scalar_t scalar_t;
	typedef typename traits_t::transpose_t transpose_t;
	typedef typename traits_t::column_t column_t;
	typedef typename traits_t::row_t row_t;
	static constexpr int rows = traits_t::rows;
	static constexpr int cols = traits_t::cols;
};

template <typename Scalar, int Rows, int Cols>
struct traits<Matrix<Scalar, Rows, Cols>> {
	typedef Matrix<Scalar, Rows, Cols> matrix_t;
	typedef TransposeMatrix<matrix_t> transpose_t;
	typedef Column<matrix_t> column_t;
	typedef Row<matrix_t> row_t;
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

    typedef typename traits<Derived>::row_t column_t;
    typedef typename traits<Derived>::column_t row_t;
};

template <class Derived>
struct traits<Column<Derived>> {
	typedef Column<Derived> self_t;
	typedef typename traits<Derived>::scalar_t scalar_t;
	typedef typename traits<Derived>::transpose_t transpose_t;
	typedef self_t column_t;
    typedef Row<self_t> row_t;

	static constexpr int rows = traits<Derived>::rows;
	static constexpr int cols = 1;
};

template <class Derived>
struct traits<Row<Derived>> {
	typedef Row<Derived> self_t;
	typedef typename traits<Derived>::scalar_t scalar_t;
	typedef typename traits<Derived>::transpose_t transpose_t;
	typedef self_t row_t;
    typedef Column<self_t> column_t;

	static constexpr int rows = 1;
	static constexpr int cols = traits<Derived>::cols;
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
	typedef Column<matrix_t> column_t;
	typedef Row<matrix_t> row_t;
};

template <class Derived, class Scalar, class BinaryOp>
struct traits<CmptWiseMatrixScalarOpResult<Derived, Scalar, BinaryOp>> {
	typedef typename traits<Derived>::scalar_t scalar_t;
	typedef CmptWiseMatrixScalarOpResult<Derived, Scalar, BinaryOp> matrix_t;

	static constexpr int rows = traits<Derived>::rows;
	static constexpr int cols = traits<Derived>::cols;

	typedef TransposeMatrix<matrix_t> transpose_t;
	typedef Column<matrix_t> column_t;
	typedef Row<matrix_t> row_t;
};

template <class Derived, class OtherDerived>
struct traits<MatrixProductResult<Derived, OtherDerived>> {
	typedef typename traits<Derived>::scalar_t scalar_t;
	typedef MatrixProductResult<Derived, OtherDerived> matrix_t;

	static constexpr int rows = traits<Derived>::rows;
	static constexpr int cols = traits<OtherDerived>::cols;

	typedef TransposeMatrix<matrix_t> transpose_t;
	typedef Column<matrix_t> column_t;
	typedef Row<matrix_t> row_t;
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

template <class Derived, class OtherDerived>
void assert_same_shape(const MatrixBase<Derived> &a, const MatrixBase<OtherDerived> &b) {
	static_assert(traits<Derived>::rows == traits<OtherDerived>::rows,
				  "Matrix dimension mismatch (rows)");
	static_assert(traits<Derived>::cols == traits<OtherDerived>::cols,
				  "Matrix dimension mismatch (cols)");
}

template <class Derived, class OtherDerived>
void assert_same_scalar_type(const MatrixBase<Derived> &a, const MatrixBase<OtherDerived> &b) {
	static_assert(std::is_same<typename traits<Derived>::scalar_t, typename traits<OtherDerived>::scalar_t>::value,
			      "Matrix scalar type mismatch");
}

template <class Derived, class OtherDerived>
void assert_multipliable(const MatrixBase<Derived> &a, const MatrixBase<OtherDerived> &b) {
	static_assert(traits<Derived>::cols == traits<OtherDerived>::rows,
				  "Matrix mutliplication dimension mismatch");
}

template <class Derived>
void assert_square(const MatrixBase<Derived> &a) {
	static_assert(traits<Derived>::rows == traits<Derived>::cols,
				  "Non-square matrix");
}

template <int Cols, class Derived>
void assert_num_cols(const MatrixBase<Derived> &a) {
    static_assert(traits<Derived>::cols == Cols, "Wrong num cols");
}


template <int Rows, class Derived>
void assert_num_rows(const MatrixBase<Derived> &a) {
    static_assert(traits<Derived>::cols == Rows, "Wrong num rows");
}

template <int Dim, class Derived>
void assert_is_vector(const MatrixBase<Derived> &a) {
    assert_num_cols<1>(a);
    assert_num_rows<Dim>(a);
}

template <class Derived>
class MatrixBase {
 public:
 	Derived &derived() { return static_cast<Derived &>(*this); }
	const Derived &derived() const { return static_cast<const Derived &>(*this); }

	typedef typename traits<Derived>::scalar_t scalar_t;
	typedef typename traits<Derived>::transpose_t transpose_t;
	typedef typename traits<Derived>::column_t column_t;
	typedef typename traits<Derived>::row_t row_t;

	int rows() const { return derived().rows(); }
	int cols() const { return derived().cols(); }

	scalar_t operator()(int i, int j) const { return derived()(i, j); }
	scalar_t &operator()(int i, int j) { return derived()(i, j); }

	transpose_t transpose() const { return derived().transpose(); }

	const column_t col(int j) const { return derived().col(j); }
	column_t col(int j) { return derived().col(j); }

	const row_t row(int i) const { return derived().row(i); }
	row_t row(int i) { return derived().row(i); }

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
        assert_square(a);
        for (int i = 0; i < a.rows(); ++i) {
            for (int j = 0; j < a.cols(); ++j) {
                a(i, j) = static_cast<int>(i == j);
            }
        }
        return a;
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
		assert_same_shape(*this, other);
		this->assign(other);
	}
};

template <typename Scalar, int Rows, int Cols>
class Matrix : public MatrixBase<Matrix<Scalar, Rows, Cols>> {
 public:
 	typedef Matrix<Scalar, Rows, Cols> self_t;
 	typedef typename traits<self_t>::scalar_t scalar_t;
 	typedef typename traits<self_t>::transpose_t transpose_t;
 	typedef typename traits<self_t>::column_t column_t;
 	typedef typename traits<self_t>::row_t row_t;

	Matrix() = default;

	scalar_t operator()(int i, int j) const {
        assert_in_bounds(i, j, *this);
		return data_[Cols * i + j];
	}

	scalar_t &operator()(int i, int j) {
        assert_in_bounds(i, j, *this);
		return data_[Cols * i + j];
	}

    template <typename = typename std::enable_if<Cols == 1>>
    scalar_t &operator()(int i) {
        assert_in_bounds(i, 0, *this);
        return data_[Cols * i + 0];
    }

    template <typename = typename std::enable_if<Cols == 1>>
    scalar_t operator()(int i) const {
        assert_in_bounds(i, 0, *this);
        return data_[Cols * i + 0];
    }

	transpose_t transpose() const { return transpose_t(*this); }

    int rows() const { return Rows; }
    int cols() const { return Cols; }

	const column_t col(int j) const { return column_t(*this, j); }
	column_t col(int j) { return column_t(*this, j); }

	const row_t row(int i) const { return row_t(*this, i); }
	row_t row(int i) { return row_t(*this, i); }

	template <class OtherDerived>
	Matrix(const MatrixBase<OtherDerived> &other) {
		this->safe_assign(other);
	}
 private:
 	Scalar data_[Rows * Cols];
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

	int rows() const { return mat_.cols(); }
	int cols() const { return mat_.rows(); }

	scalar_t operator()(int i, int j) const { return mat_(j, i); }

	transpose_t transpose() const { return mat_; }
 private:
 	const Derived &mat_;
};

template <class Derived, class OtherDerived, class BinaryOp>
class CmptWiseBinaryMatrixOpResult : public MatrixBase<CmptWiseBinaryMatrixOpResult<Derived, OtherDerived, BinaryOp>> {
 public:
 	CmptWiseBinaryMatrixOpResult(const Derived &a, const OtherDerived &b) : a_(a), b_(b) {
		assert_same_shape(a, b);
		assert_same_scalar_type(a, b);
	}

	const BinaryOp op = BinaryOp();

	typedef CmptWiseBinaryMatrixOpResult<Derived, OtherDerived, BinaryOp> self_t;
	typedef typename traits<self_t>::scalar_t scalar_t;
	typedef typename traits<self_t>::transpose_t transpose_t;

	int rows() const { return a_.rows(); }
	int cols() const { return a_.cols(); }

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

	int rows() const { return a_.rows(); }
	int cols() const { return a_.cols(); }

	scalar_t operator()(int i, int j) const { return op(a_(i, j), static_cast<scalar_t>(s_)); }
 private:
	const Derived &a_;
 	const Scalar s_;
};

template <class Derived, class OtherDerived>
class MatrixProductResult : public MatrixBase<MatrixProductResult<Derived, OtherDerived>> {
 public:
 	MatrixProductResult(const Derived &a, const OtherDerived &b) : a_(a), b_(b) {
        std::cout << a_ << std::endl;
		assert_same_scalar_type(a, b);
		assert_multipliable(a, b);
	}

	typedef MatrixProductResult<Derived, OtherDerived> self_t;
	typedef typename traits<self_t>::scalar_t scalar_t;

	int rows() const { return a_.rows(); }
	int cols() const { return b_.cols(); }

	scalar_t operator()(int i, int j) const {
		scalar_t sum = 0;
        std::cout << a_ << std::endl;
		for (int k = 0; k < a_.cols(); ++k) {
			sum += a_(i, k) * b_(k, j);
		}
		return sum;
	}
 private:
	const Derived &a_;
 	const OtherDerived &b_;
};

template <class Derived>
class Column : public MatrixBase<Column<Derived>> {
 public:
     Column(const Derived &a, int j) : const_a_(a), j_(j) {}
     Column(Derived &a, int j) : const_a_(a), a_(&a), j_(j) {}

     typedef Column<Derived> self_t;
     typedef typename traits<self_t>::scalar_t scalar_t;

     int rows() const { return const_a_.rows(); }
     int cols() const { return 1; }

     scalar_t operator()(int i, int j) const {
         assert(j == 0);
         return const_a_(i, j_);
     }

     scalar_t &operator()(int i, int j) {
         assert(j == 0);
         return (*a_)(i, j_);
     }

     scalar_t operator()(int i) const {
         return const_a_(i, j_);
     }

     scalar_t &operator()(int i) {
         return (*a_)(i, j_);
     }
 private:
     const Derived &const_a_;
     Derived *a_;
     int j_;
};

template <typename Derived>
class Row : public MatrixBase<Row<Derived>> {
 public:
     Row(const Derived &a, int i) : const_a_(a), i_(i) {}
     Row(Derived &a, int i) : const_a_(a), a_(&a), i_(i) {}

     typedef Row<Derived> self_t;
     typedef typename traits<self_t>::scalar_t scalar_t;

     int cols() const { return const_a_.cols(); }
     int rows() const { return 1; }

     scalar_t operator()(int i, int j) const {
         assert(i == 0);
         return const_a_(i_, j);
     }

     scalar_t &operator()(int i, int j) {
         assert(i == 0);
         return (*a_)(i_, j);
     }

     scalar_t operator()(int j) const {
         return const_a_(i_, j);
     }

     scalar_t &operator()(int j) {
         return (*a_)(i_, j);
     }
 private:
     const Derived &const_a_;
     Derived *a_;
     int i_;
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

}  // namspace sfm
