#include "linalg/matrix.hh"

#include <cassert>
#include <iostream>

namespace sfm {
namespace linalg {

template <int ROWS, int COLS>
Matrix<double, ROWS, COLS> make_matrix() {
    Matrix<double, ROWS, COLS> m;
    for (int i = 0; i < ROWS; ++i) {
        for (int j = 0; j < COLS; ++j) {
            m(i, j) = static_cast<double>(ROWS  * i - COLS * j);
        }
    }
    return m;
}

void rows_cols_test() {
    constexpr int ROWS = 10;
    constexpr int COLS = 5;
    const auto m = make_matrix<ROWS, COLS>();
    assert(m.rows() == ROWS);
    assert(m.cols() == COLS);
}

void fill_test() {
    constexpr int ROWS = 4;
    constexpr int COLS = 3;
    Matrix<double, ROWS, COLS> all_7s;

    constexpr double FILL_VALUE = 7.0;
    all_7s.fill(FILL_VALUE);
    for (int i = 0; i < ROWS; ++i) {
        for (int j = 0; j < COLS; ++j) {
            assert(all_7s(i, j) == FILL_VALUE);
        }
    }
}

void zero_matrix_test() {
    constexpr int ROWS = 3;
    constexpr int COLS = 5;
    const Matrix<double, ROWS, COLS> zm =
        Matrix<double, ROWS, COLS>::zero();
    for (int i = 0; i < ROWS; ++i) {
        for (int j = 0; j < COLS; ++j) {
            assert(zm(i, j) == 0.0);
        }
    }
}

void identity_matrix_test() {
    constexpr int DIM = 9;
    const Matrix<double, DIM, DIM> id =
        Matrix<double, DIM, DIM>::identity();
    for (int i = 0; i < DIM; ++i) {
        for (int j = 0; j < DIM; ++j) {
            if (i == j) {
                assert(id(i, j) == 1.0);
            } else {
                assert(id(i, j) == 0.0);
            }
        }
    }
}

void transpose_test() {
    constexpr int ROWS = 11;
    constexpr int COLS = 8;
    const auto m = make_matrix<ROWS, COLS>();
    const Matrix<double, COLS, ROWS> mt = m.transpose();
    assert(m.rows() == mt.cols());
    assert(m.cols() == mt.rows());
    for (int i = 0; i < ROWS; ++i) {
        for (int j = 0; j < COLS; ++j) {
            assert(m(i, j) == mt(j, i));
        }
    }

    const Matrix<double, ROWS, COLS> mtt = m.transpose().transpose();
    for (int i = 0; i < ROWS; ++i) {
        for (int j = 0; j < COLS; ++j) {
            assert(m(i, j) == mtt(i, j));
        }
    }
}

void col_test() {
    constexpr int ROWS = 9;
    constexpr int COLS = 3;
    auto m = make_matrix<ROWS, COLS>();

    constexpr int COL = 1;
    const Matrix<double, ROWS, 1> col = m.col(COL);

    for (int i = 0; i < ROWS; ++i) {
        assert(col(i) == m(i, COL));
    }

    const auto to_add = make_matrix<ROWS, 1>();

    m.col(COL) = m.col(COL);
    m.col(COL) += to_add;
    for (int i = 0; i < ROWS; ++i) {
        assert(m(i, COL) == col(i) + to_add(i, 0));
    }
}

void row_test() {
    constexpr int ROWS = 3;
    constexpr int COLS = 7;
    auto m = make_matrix<ROWS, COLS>();

    constexpr int ROW = 2;
    const Matrix<double, 1, COLS> row = m.row(ROW);

    for (int j = 0; j < COLS; ++j) {
        assert(row(j) == m(ROW, j));
    }

    const auto to_add = make_matrix<1, COLS>();

    m.row(ROW) = m.row(ROW);
    m.row(ROW) += to_add;
    for (int j = 0; j < COLS; ++j) {
        assert(m(ROW, j) == row(j) + to_add(j));
    }
}

void cmptwise_binary_op_test() {
    constexpr int ROWS = 5;
    constexpr int COLS = 3;
    const auto m = make_matrix<ROWS, COLS>();
    const auto n = make_matrix<ROWS, COLS>();
    const Matrix<double, ROWS, COLS> sum = m + n;
    const Matrix<double, ROWS, COLS> diff = m - n;
    for (int i = 0; i < ROWS; ++i) {
        for (int j = 0; j < COLS; ++j) {
            assert(sum(i, j) == m(i, j) + n(i, j));
            assert(diff(i, j) == m(i, j) - n(i, j));
        }
    }
}

void cmptwise_scalar_op_test() {
    constexpr int ROWS = 7;
    constexpr int COLS = 11;

    const auto m = make_matrix<ROWS, COLS>();

    constexpr double S = 2.0;
    const Matrix<double, ROWS, COLS> m_times_s = m * S;
    const Matrix<double, ROWS, COLS> m_by_s = m / S;

    for (int i = 0; i < ROWS; ++i) {
        for (int j = 0; j < COLS; ++j) {
            assert(m_times_s(i, j) == m(i, j) * S);
            assert(m_by_s(i, j) == m(i, j) / S);
        }
    }
}

void matrix_product_test() {
    Matrix<int, 2, 3> m;
    m(0, 0) = 1; m(0, 1) = 2; m(0, 2) = 3;
    m(1, 0) = 4; m(1, 1) = 5; m(1, 2) = 6;

    Matrix<int, 3, 2> n;
    n(0, 0) = 7;  n(0, 1) = 8;
    n(1, 0) = 9;  n(1, 1) = 10;
    n(2, 0) = 11; n(2, 1) = 12;

    const Matrix<int, 2, 2> mn = m * n;
    assert(mn(0, 0) == 58);  assert(mn(0, 1) == 64);
    assert(mn(1, 0) == 139); assert(mn(1, 1) == 154);
}

}
}

int main() {
    std::cout << "rows_cols_test..." << std::endl;
    sfm::linalg::rows_cols_test();

    std::cout << "zero_matrix_test..." << std::endl;
    sfm::linalg::zero_matrix_test();

    std::cout << "identity_matrix_test..." << std::endl;
    sfm::linalg::identity_matrix_test();

    std::cout << "transpose_test..." << std::endl;
    sfm::linalg::transpose_test();

    std::cout << "col_test..." << std::endl;
    sfm::linalg::col_test();

    std::cout << "row_test..." << std::endl;
    sfm::linalg::row_test();

    std::cout << "cmptwise_scalar_op_test..." << std::endl;
    sfm::linalg::cmptwise_scalar_op_test();

    std::cout << "cmptwise_binary_op_test..." << std::endl;
    sfm::linalg::cmptwise_binary_op_test();

    std::cout << "matrix_product_test..." << std::endl;
    sfm::linalg::matrix_product_test();

    std::cout << "All tests passed!" << std::endl;
    return 0;
}
