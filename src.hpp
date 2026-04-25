#ifndef SRC_HPP
#define SRC_HPP

#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <utility>
#include <cstdlib>
#include <cmath>
#include <exception>

#include "fraction.hpp"

class matrix {
private:

    int m, n;
    fraction **data;

    void allocate(int m_, int n_) {
        m = m_;
        n = n_;
        if (m == 0 || n == 0) {
            data = nullptr;
            return;
        }
        data = new fraction *[m];
        for (int i = 0; i < m; ++i) {
            data[i] = new fraction[n];
            for (int j = 0; j < n; ++j) data[i][j] = fraction(0);
        }
    }

public:

    matrix() {
        m = n = 0;
        data = nullptr;
    }

    matrix(int m_, int n_) { allocate(m_, n_); }

    matrix(const matrix &obj) {
        allocate(obj.m, obj.n);
        if (data) {
            for (int i = 0; i < m; ++i)
                for (int j = 0; j < n; ++j)
                    data[i][j] = obj.data[i][j];
        }
    }

    matrix(matrix &&obj) noexcept {
        m = obj.m; n = obj.n; data = obj.data;
        obj.m = obj.n = 0; obj.data = nullptr;
    }

    ~matrix() {
        if (data) {
            for (int i = 0; i < m; ++i) delete [] data[i];
            delete [] data;
        }
        data = nullptr; m = n = 0;
    }

    matrix &operator=(const matrix &obj) {
        if (this == &obj) return *this;
        // cleanup current
        if (data) {
            for (int i = 0; i < m; ++i) delete [] data[i];
            delete [] data;
        }
        allocate(obj.m, obj.n);
        if (data) {
            for (int i = 0; i < m; ++i)
                for (int j = 0; j < n; ++j)
                    data[i][j] = obj.data[i][j];
        }
        return *this;
    }

    // i: 1-based row, j: 0-based col
    fraction &operator()(int i, int j) {
        if (i < 1 || i > m || j < 0 || j >= n || data == nullptr) {
            throw matrix_error();
        }
        return data[i - 1][j];
    }

    friend matrix operator*(const matrix &lhs, const matrix &rhs) {
        if (lhs.n != rhs.m || lhs.m == 0 || lhs.n == 0 || rhs.n == 0) {
            throw matrix_error();
        }
        matrix res(lhs.m, rhs.n);
        for (int i = 0; i < lhs.m; ++i) {
            for (int k = 0; k < lhs.n; ++k) {
                const fraction &lik = lhs.data[i][k];
                if (lik == fraction(0)) continue;
                for (int j = 0; j < rhs.n; ++j) {
                    res.data[i][j] = res.data[i][j] + lik * rhs.data[k][j];
                }
            }
        }
        return res;
    }

    matrix transposition() {
        if (m == 0 || n == 0 || data == nullptr) throw matrix_error();
        matrix t(n, m);
        for (int i = 0; i < m; ++i)
            for (int j = 0; j < n; ++j)
                t.data[j][i] = data[i][j];
        return t;
    }

    fraction determination() {
        if (m == 0 || n == 0 || m != n || data == nullptr) throw matrix_error();
        int sz = n;
        // make a copy
        matrix a(*this);
        bool neg = false;
        fraction det(1);
        for (int col = 0; col < sz; ++col) {
            int pivot = col;
            while (pivot < sz && a.data[pivot][col] == fraction(0)) pivot++;
            if (pivot == sz) return fraction(0);
            if (pivot != col) {
                std::swap(a.data[pivot], a.data[col]);
                neg = !neg;
            }
            fraction piv = a.data[col][col];
            for (int row = col + 1; row < sz; ++row) {
                if (a.data[row][col] == fraction(0)) continue;
                fraction factor = a.data[row][col] / piv;
                for (int c = col; c < sz; ++c) {
                    a.data[row][c] = a.data[row][c] - factor * a.data[col][c];
                }
            }
        }
        for (int i = 0; i < sz; ++i) det = det * a.data[i][i];
        if (neg) det = fraction(0) - det; // multiply by -1
        return det;
    }

    // helpers
    int rows() const { return m; }
    int cols() const { return n; }
};

class resistive_network {
private:

    int interface_size, connection_size;
    matrix adjacency, conduction;

    std::vector<int> edge_from;
    std::vector<int> edge_to;
    std::vector<fraction> edge_res;
    std::vector<fraction> edge_cond; // conductance 1/r

    static std::vector<fraction> solve_linear(matrix A, const std::vector<fraction> &b) {
        int n = A.rows();
        if (n == 0 || A.cols() != n || (int)b.size() != n) throw matrix_error();
        std::vector<fraction> rhs = b;
        for (int col = 0; col < n; ++col) {
            // find pivot row with nonzero in this column
            int piv_row = -1;
            for (int r = col; r < n; ++r) {
                if (!(A(r + 1, col) == fraction(0))) { piv_row = r; break; }
            }
            if (piv_row == -1) throw matrix_error();
            if (piv_row != col) {
                for (int j = 0; j < n; ++j) std::swap(A(col + 1, j), A(piv_row + 1, j));
                std::swap(rhs[col], rhs[piv_row]);
            }
            fraction piv = A(col + 1, col);
            for (int r = col + 1; r < n; ++r) {
                if (A(r + 1, col) == fraction(0)) continue;
                fraction factor = A(r + 1, col) / piv;
                for (int c = col; c < n; ++c) {
                    A(r + 1, c) = A(r + 1, c) - factor * A(col + 1, c);
                }
                rhs[r] = rhs[r] - factor * rhs[col];
            }
        }
        std::vector<fraction> x(n, fraction(0));
        for (int i = n - 1; i >= 0; --i) {
            fraction sum(0);
            for (int j = i + 1; j < n; ++j) sum = sum + A(i + 1, j) * x[j];
            x[i] = (rhs[i] - sum) / A(i + 1, i);
        }
        return x;
    }

public:

    resistive_network(int interface_size_, int connection_size_, int from[], int to[], fraction resistance[])
            : interface_size(interface_size_), connection_size(connection_size_),
              adjacency(connection_size_, interface_size_), conduction(connection_size_, connection_size_) {
        edge_from.resize(connection_size_);
        edge_to.resize(connection_size_);
        edge_res.resize(connection_size_);
        edge_cond.resize(connection_size_);
        for (int k = 0; k < connection_size_; ++k) {
            int u = from[k];
            int v = to[k];
            edge_from[k] = u;
            edge_to[k] = v;
            edge_res[k] = resistance[k];
            edge_cond[k] = fraction(1) / resistance[k];
            // set incidence row k: +1 at u, -1 at v
            adjacency(k + 1, u - 1) = adjacency(k + 1, u - 1) + fraction(1);
            adjacency(k + 1, v - 1) = adjacency(k + 1, v - 1) - fraction(1);
            // conduction diagonal
            conduction(k + 1, k) = conduction(k + 1, k) + edge_cond[k];
        }
    }

    ~resistive_network() = default;

    fraction get_equivalent_resistance(int interface_id1, int interface_id2) {
        // L = A^T C A
        matrix L = adjacency.transposition() * conduction * adjacency;
        int n = interface_size;
        // build reduced matrix by removing last node
        matrix Lr(n - 1, n - 1);
        for (int i = 0; i < n - 1; ++i)
            for (int j = 0; j < n - 1; ++j)
                Lr(i + 1, j) = L(i + 1, j);
        std::vector<fraction> b(n - 1, fraction(0));
        if (interface_id1 < n) b[interface_id1 - 1] = b[interface_id1 - 1] + fraction(1);
        if (interface_id2 < n) b[interface_id2 - 1] = b[interface_id2 - 1] - fraction(1);
        auto U = solve_linear(Lr, b);
        fraction ui = (interface_id1 < n) ? U[interface_id1 - 1] : fraction(0);
        fraction uj = (interface_id2 < n) ? U[interface_id2 - 1] : fraction(0);
        return ui - uj;
    }

    fraction get_voltage(int id, fraction current[]) {
        matrix L = adjacency.transposition() * conduction * adjacency;
        int n = interface_size;
        matrix Lr(n - 1, n - 1);
        for (int i = 0; i < n - 1; ++i)
            for (int j = 0; j < n - 1; ++j)
                Lr(i + 1, j) = L(i + 1, j);
        std::vector<fraction> b(n - 1, fraction(0));
        for (int i = 0; i < n - 1; ++i) b[i] = current[i];
        auto U = solve_linear(Lr, b);
        return U[id - 1];
    }

    fraction get_power(fraction voltage[]) {
        fraction total(0);
        for (int k = 0; k < connection_size; ++k) {
            int u = edge_from[k];
            int v = edge_to[k];
            fraction du = voltage[u - 1] - voltage[v - 1];
            total = total + edge_cond[k] * du * du;
        }
        return total;
    }
};


#endif //SRC_HPP
