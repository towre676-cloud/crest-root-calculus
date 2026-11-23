// Exact periodic cyclic Jacobi characteristic polynomial, GMP integers.
// Input: length, iteration count, diagonal entries, cyclic bond entries.
// Standard transfer matrices and Chebyshev recurrence; see proof.md.
#include <gmpxx.h>
#include <algorithm>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

using Z = mpz_class;
using Poly = std::vector<Z>;

static void trim(Poly& p) {
    while (p.size() > 1 && p.back() == 0) p.pop_back();
}
static Poly linear(const Poly& p, const Z& a) {
    Poly out(p.size() + 1);
    for (size_t j = 0; j < p.size(); ++j) {
        out[j + 1] += p[j];
        mpz_submul(out[j].get_mpz_t(), p[j].get_mpz_t(), a.get_mpz_t());
    }
    trim(out); return out;
}
static void add_scaled(Poly& p, const Poly& q, const Z& c) {
    p.resize(std::max(p.size(), q.size()));
    for (size_t j = 0; j < q.size(); ++j)
        mpz_addmul(p[j].get_mpz_t(), q[j].get_mpz_t(), c.get_mpz_t());
    trim(p);
}
static Poly scaled(Poly p, const Z& c) {
    for (auto& v : p) v *= c;
    trim(p); return p;
}
static Poly multiply(const Poly& p, const Poly& q) {
    Poly out(p.size() + q.size() - 1);
    for (size_t i = 0; i < p.size(); ++i) if (p[i] != 0)
        for (size_t j = 0; j < q.size(); ++j) if (q[j] != 0)
            mpz_addmul(out[i + j].get_mpz_t(), p[i].get_mpz_t(), q[j].get_mpz_t());
    trim(out); return out;
}
static size_t period(const std::vector<Z>& a, const std::vector<Z>& b) {
    std::vector<size_t> prefix(a.size());
    for (size_t i = 1; i < a.size(); ++i) {
        size_t k = prefix[i - 1];
        while (k && (a[i] != a[k] || b[i] != b[k])) k = prefix[k - 1];
        if (a[i] == a[k] && b[i] == b[k]) ++k;
        prefix[i] = k;
    }
    size_t m = a.size() - prefix.back();
    return a.size() % m ? a.size() : m;
}
static Poly solve(const std::vector<Z>& a, const std::vector<Z>& b, size_t& m) {
    if (a.size() < 3 || a.size() != b.size()) throw std::runtime_error("bad length");
    for (const auto& v : b) if (v == 0) throw std::runtime_error("zero bond");
    m = period(a,b);
    const size_t n = a.size() / m;
    Poly p{1}, q{0}, r{0}, s{1};
    Z product = 1;
    for (size_t i = 0; i < m; ++i) {
        Poly newp = linear(p,a[i]), newq = linear(q,a[i]);
        Z c = -b[(i + m - 1) % m];
        add_scaled(newp,r,c); add_scaled(newq,s,c);
        r = scaled(std::move(p), b[i]); s = scaled(std::move(q), b[i]);
        p = std::move(newp); q = std::move(newq); product *= b[i];
    }
    add_scaled(p,s,Z(1));
    const Poly d = p;
    Poly previous{2}, current = d;
    Z minus_square = -product*product;
    for (size_t k = 2; k <= n; ++k) {
        Poly next = multiply(d,current);
        add_scaled(next,previous,minus_square);
        previous = std::move(current); current = std::move(next);
    }
    Z bn; mpz_pow_ui(bn.get_mpz_t(), product.get_mpz_t(), n);
    current[0] -= 2*bn;
    return current;
}
static double cpu_time() {
    struct timespec t;
    if (clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&t)) throw std::runtime_error("clock failed");
    return t.tv_sec + t.tv_nsec*1e-9;
}
int main() {
    try {
        size_t length, iterations;
        if (!(std::cin >> length >> iterations) || length < 3 || iterations < 1)
            throw std::runtime_error("bad header");
        std::vector<Z> a(length), b(length);
        for (auto& v : a) if (!(std::cin >> v)) throw std::runtime_error("bad diagonal");
        for (auto& v : b) if (!(std::cin >> v)) throw std::runtime_error("bad bond");
        Poly result; size_t m=0;
        const double start = cpu_time();
        for (size_t k = 0; k < iterations; ++k) result=solve(a,b,m);
        const double seconds = cpu_time()-start;
        std::cout << "{\"period\":" << m << ",\"iterations\":" << iterations
                  << ",\"cpu_seconds\":" << std::setprecision(12) << seconds
                  << ",\"coefficients\":[";
        for (size_t j=0;j<result.size();++j) {
            if (j) std::cout << ',';
            std::cout << '\"' << result[j] << '\"';
        }
        std::cout << "]}\n";
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n'; return 1;
    }
}
