#pragma once

#include "linalg/matrix.hh"

#include <ostream>

namespace sfm {
namespace linalg {

template <typename Derived>
std::ostream &operator<<(std::ostream &out, const MatrixBase<Derived> &A) {
	out << "[";
	for (int i = 0; i < A.rows(); ++i) {
		out << ((i == 0) ? "[" : " [");
		for (int j = 0; j < A.cols(); ++j) {
			out << A(i, j) << (j == A.cols() - 1 ? "" : ", ");
		}
		out << ((i == A.rows() - 1) ? "]" : "]\n");
	}
    out << "]\n";
	return out;
}

}
}
