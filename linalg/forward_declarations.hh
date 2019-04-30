#pragma once

namespace sfm {

template <typename T>
struct traits;

template <class Derived>
class MatrixBase;

template <typename Scalar, int Rows, int Cols>
class Matrix;

template <class Mat>
class Column;

template <class Mat>
class Row;

template <class Derived>
class TransposeMatrix;

template <class Derived, class OtherDerived, class BinaryOp>
class CmptWiseBinaryMatrixOpResult;

template <class Derived, class Scalar, class BinaryOp>
class CmptWiseMatrixScalarOpResult;

template <class Derived, class OtherDerived>
class MatrixProductResult;

template <class Derived>
class Column;

template <class Derived>
class Row;

}
