#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <vector>

using i64 = long long;
using Point = std::array<int, 3>;

i64 isqrt(i64 n) {
    auto r = static_cast<i64>(std::sqrt(static_cast<double>(n)));
    while ((r+1)*(r+1) <= n) ++r;
    while (r*r > n) --r;
    return r;
}

i64 floor_div(i64 n, i64 d) {
    return n/d - (n%d < 0);
}

int mod3(i64 n) {
    return static_cast<int>((n%3+3)%3);
}

Point normal(Point p, Point q) {
    return {p[1]*q[2]-p[2]*q[1], p[2]*q[0]-p[0]*q[2], p[0]*q[1]-p[1]*q[0]};
}

struct Parent {
    i64 det, a, b, c, u, v, w;
};

struct Counts {
    int first = 0;
    int excess = 0;
    int witness = 0;
    bool product = false;
};

template<class Visit>
void visit_squares(Parent g, const std::vector<int>& roots, int bound, Visit visit) {
    const auto [det,a,b,c,u,v,w] = g;
    const i64 minor = a*b-u*u;
    const i64 coupling = a*w-u*v;
    const i64 transverse = a*c-v*v;
    const i64 cutoff = static_cast<i64>(bound)*bound;
    const i64 zmax = isqrt(minor*cutoff/det);
    for (i64 z = 0; z <= zmax; ++z) {
        const i64 yzcenter = coupling*z;
        const i64 yradius = isqrt(a*(minor*cutoff-det*z*z));
        const i64 ymin = z == 0 ? 0 : -floor_div(yradius+yzcenter,minor);
        const i64 ymax = floor_div(yradius-yzcenter,minor);
        for (i64 y = ymin; y <= ymax; ++y) {
            const i64 xcenter = u*y+v*z;
            const i64 minimum = minor*y*y+2*yzcenter*y+transverse*z*z;
            const i64 remaining = a*cutoff-minimum;
            if (remaining < 0) continue;
            const i64 xradius = isqrt(remaining);
            const i64 xmin = y == 0 && z == 0 ? 1 : -floor_div(xradius+xcenter,a);
            const i64 xmax = floor_div(xradius-xcenter,a);
            const int residue_base = 3*mod3(y)+9*mod3(z);
            for (i64 x = xmin; x <= xmax; ++x) {
                const i64 norm = (a*x+2*xcenter)*x+b*y*y+2*w*y*z+c*z*z;
                const int m = roots[norm];
                if (m && !visit(m,mod3(x)+residue_base)) return;
            }
        }
    }
}

Counts compare(Parent g, Point l, Point r, const std::vector<int>& roots, int bound) {
    std::array<int, 27> sign;
    for (int z = 0; z < 3; ++z)
        for (int y = 0; y < 3; ++y)
            for (int x = 0; x < 3; ++x) {
                const bool left = mod3(l[0]*x+l[1]*y+l[2]*z) == 0;
                const bool right = mod3(r[0]*x+r[1]*y+r[2]*z) == 0;
                sign[x+3*y+9*z] = left == right ? 2*left : static_cast<int>(left)-right;
            }
    Counts result;
    std::vector<int> differences(bound+1);
    visit_squares(g,roots,bound,[&](int m,int residue) {
        const int membership = sign[residue];
        if (membership == 2) {
            if (m == 1) {
                result.product = true;
                return false;
            }
            if ((!result.witness || m < result.witness) && std::gcd<i64>(m, 6*g.det) == 1)
                result.witness = m;
        } else if (membership) {
            if (m%3 != 0) std::abort();
            differences[m] += 2*membership;
        }
        return true;
    });
    for (int m = 3; m <= bound; m += 3)
        if (differences[m]) {
            result.first = m;
            result.excess = differences[m];
            break;
        }
    return result;
}

int main(int argc, char** argv) {
    if ((argc != 2 && argc != 3) || std::atoi(argv[1]) < 3 || std::atoi(argv[1]) > 8192) {
        std::cerr << "usage: shells MAX_DEGREE (3..8192) [--theta]\n";
        return 2;
    }
    const int cap = std::atoi(argv[1]);
    std::vector<int> roots(cap*cap+1);
    for (int m = 1; m <= cap; ++m) roots[m*m] = m;
    Parent g;
    using Edge = std::array<Point, 2>;
    while (std::cin >> g.det >> g.a >> g.b >> g.c >> g.u >> g.v >> g.w) {
        if (argc == 3) {
            int bound;
            if (!(std::cin >> bound) || bound < 1 || bound > cap) return 2;
            std::vector<int> counts(bound+1);
            visit_squares(g,roots,bound,[&](int m,int) {
                counts[m] += 2;
                return true;
            });
            for (int m = 1; m <= bound; ++m) std::cout << counts[m] << ' ';
            std::cout << '\n';
            continue;
        }
        std::array<Edge, 2> edges;
        for (auto& edge : edges)
            for (auto& p : edge)
                for (int& x : p)
                    if (!(std::cin >> x)) return 2;
        Point l = normal(edges[0][0],edges[0][1]);
        Point r = normal(edges[1][0],edges[1][1]);
        int bound = std::min<int>(cap, 6+2*isqrt(g.c));
        Counts result;
        for (;;) {
            result = compare(g,l,r,roots,bound);
            if ((result.first && result.witness) || result.product || bound == cap) break;
            bound = std::min(cap, 2*bound);
        }
        std::cout << result.first << ' ' << result.excess << ' ' << result.witness << ' '
                  << result.product << ' ' << bound << '\n';
    }
}
