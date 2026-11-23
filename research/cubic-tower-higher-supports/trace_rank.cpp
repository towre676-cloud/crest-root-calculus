// Exact prime-field fibre counts and modular row reduction for the signed tower.
// Input: one simultaneous good prime per line, checked by the Python driver.
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

constexpr int stages = 10;
constexpr int subsets = 1 << stages;
constexpr int columns = 2 * subsets - 1;
constexpr int modulus = 127;
using Row = std::array<int, columns>;

static int residue(std::int64_t x, int p) {
    int a = x % p;
    return a < 0 ? a + p : a;
}

static Row traces(int p) {
    std::vector<int> root(p, -1), inverse(p), inverse_cube(p);
    for (int x = 0; x < p; ++x) root[std::int64_t(x)*x % p] = x;
    int nonsquare = 2;
    while (root[nonsquare] >= 0) ++nonsquare;
    inverse[1] = 1;
    for (int x = 2; x < p; ++x) inverse[x] = p - std::int64_t(p/x)*inverse[p%x] % p;
    std::array<std::vector<int>, 3> freq;
    for (auto& f : freq) f.assign(p, 0);
    for (int x = 0; x < p; ++x) {
        const auto cube = std::int64_t(x)*x % p*x % p;
        ++freq[0][cube];
        ++freq[1][(cube+x) % p];
        ++freq[2][(cube+std::int64_t(nonsquare)*x) % p];
        inverse_cube[x] = inverse[cube];
    }
    std::vector<int> type(p, 0), scale(p, 1);
    for (int a = 1; a < p; ++a) {
        type[a] = root[a] < 0 ? 2 : 1;
        const int a_square = type[a] == 2 ? std::int64_t(a)*inverse[nonsquare] % p : a;
        scale[a] = inverse_cube[root[a_square]];
    }
    std::array<std::array<int, subsets>, 2> total{};
    int base_points = 1;  // infinity
    for (int twist = 0; twist < 2; ++twist) {
        int inv_s = twist == 0 ? 1 : inverse[nonsquare];
        std::array<int, subsets> product{};
        for (int z = 0; z <= p; ++z) {
            std::array<int, stages> chars;
            if (z == p) {
                product[0] = 1;
                for (int i = 0; i < stages; ++i) {
                    const int a = residue(-std::int64_t(12*i+1)*inverse[12 % p], p);
                    const int chi = root[a] >= 0 ? 1 : -1;
                    chars[i] = 1 + (twist == 0 ? chi : -chi);
                }
            } else {
                const auto f = residue(((-4LL*z+1)*z+1)*z+1, p);
                const int w = root[std::int64_t(f)*inv_s % p];
                if (w < 0) continue;
                product[0] = w == 0 ? 1 : 2;
                if (twist == 0) base_points += product[0];
                const auto b = std::int64_t(w)*inv_s % p;
                for (int i = 0; i < stages; ++i) {
                    const int a = 3LL*(z+i)*inv_s % p;
                    chars[i] = freq[type[a]][b*scale[a] % p] - 1;
                }
            }
            for (int i = 0; i < stages; ++i) {
                const int n = 1 << i;
                for (int s = 0; s < n; ++s) product[n+s] = product[s]*chars[i];
            }
            for (int s = 1; s < subsets; ++s) total[twist][s] += product[s];
        }
    }
    Row result{};
    result[0] = p+1-base_points;
    for (int s = 1; s < subsets; ++s) {
        const auto a = -total[0][s], b = -total[1][s];
        if ((a+b) % 2 != 0) throw std::runtime_error("nonintegral signed trace");
        result[2*s-1] = (a+b)/2;
        result[2*s] = (a-b)/2;
    }
    return result;
}

int main(int argc, char** argv) {
    const bool raw = argc == 2 && std::string(argv[1]) == "--traces";
    std::vector<int> primes;
    for (int p; std::cin >> p;) {
        if (p < 5 || p > 100000) throw std::runtime_error("prime outside arithmetic bounds");
        primes.push_back(p);
    }
    std::vector<Row> rows(primes.size());
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 4)
#endif
    for (std::size_t i = 0; i < primes.size(); ++i) rows[i] = traces(primes[i]);
    if (raw) {
        for (const auto& row : rows) {
            for (const int v : row) std::cout << v << ' ';
            std::cout << '\n';
        }
        return 0;
    }
    std::array<int, modulus> inv{};
    for (int i = 1; i < modulus; ++i) {
        for (int j = 1; j < modulus; ++j) if (i*j % modulus == 1) inv[i] = j;
    }
    for (auto& row : rows) for (auto& x : row) x = residue(x, modulus);
    std::size_t rank = 0;
    for (int col = 0; col < columns && rank < rows.size(); ++col) {
        auto pivot = rank;
        while (pivot < rows.size() && rows[pivot][col] == 0) ++pivot;
        if (pivot == rows.size()) continue;
        std::swap(rows[rank], rows[pivot]);
        const int normalizer = inv[rows[rank][col]];
        for (int j = col; j < columns; ++j) rows[rank][j] = rows[rank][j]*normalizer % modulus;
        for (auto i = rank+1; i < rows.size(); ++i) {
            const int m = modulus-rows[i][col];
            if (m == modulus) continue;
            for (int j = col; j < columns; ++j) rows[i][j] = (rows[i][j]+m*rows[rank][j]) % modulus;
        }
        ++rank;
    }
    std::cout << rank << '\n';
}
