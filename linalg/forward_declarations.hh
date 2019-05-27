#pragma once

namespace sfm {
namespace linalg {

template <typename T>
struct traits;

template <class Derived>
class MatrixBase;

template <typename Scalar, int Rows, int Cols>
class Matrix;

template <class Derived, int Rows, int Cols>
class Block;

template <class Derived>
class TransposeMatrix;

template <class Derived, class OtherDerived, class BinaryOp>
class CmptWiseBinaryMatrixOpResult;

template <class Derived, class Scalar, class BinaryOp>
class CmptWiseMatrixScalarOpResult;

template <class Derived, class UnaryOp>
class CmptWiseUnaryMatrixOpResult;

template <class Derived, class OtherDerived>
class MatrixProductResult;

}
}
