// Exact representation counts for u^2+(2x+Ty)^2+|D|y^2, with a linear congruence mod ell.
#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

struct Congruence { int trace, omega_trace; };
struct Pair { int u, w, next; };

int root(int n) {
    int r = static_cast<int>(std::sqrt(n));
    while (1LL * (r+1) * (r+1) <= n) ++r;
    while (1LL * r * r > n) --r;
    return r;
}

int main() {
    int abs_discriminant, trace_omega, prime, form_count, index_count;
    if (!(std::cin >> abs_discriminant >> trace_omega >> prime >> form_count >> index_count)
        || abs_discriminant < 3 || (trace_omega != 0 && trace_omega != 1)
        || prime < 3 || form_count < 1 || index_count < 1) return 1;
    std::vector<Congruence> forms(form_count);
    for (auto &f : forms) if (!(std::cin >> f.trace >> f.omega_trace)) return 1;
    std::vector<int> indices(index_count);
    for (auto &n : indices) if (!(std::cin >> n) || n < 1) return 1;
    const int max = *std::max_element(indices.begin(), indices.end());
    std::vector<int> head(max+1, -1);
    std::vector<Pair> pairs;
    for (int u = 0, lim = root(max); u <= lim; ++u) {
        for (int w = 0, last = root(max-u*u); w <= last; ++w) {
            const int n = u*u+w*w;
            pairs.push_back({u,w,head[n]});
            head[n] = static_cast<int>(pairs.size())-1;
        }
    }
    for (const int n : indices) {
        if (n % 4 == 2 || n % 4 == 3) {
            for (int i=0; i<form_count; ++i) std::cout << "0 ";
            std::cout << '\n';
            continue;
        }
        std::vector<long long> counts(form_count, 0);
        for (int y=0, bound=root(n/abs_discriminant); y<=bound; ++y) {
            // Central negation pairs all y>0 solutions with y<0 solutions.
            const long long weight = y == 0 ? 1 : 2;
            for (int at=head[n-abs_discriminant*y*y]; at!=-1; at=pairs[at].next) {
                const auto [u,w,next] = pairs[at];
                if ((w-trace_omega*y)%2 != 0) continue;
                for (int su=1; su >= (u == 0 ? 1 : -1); su-=2)
                    for (int sw=1; sw >= (w == 0 ? 1 : -1); sw-=2) {
                        const int x = (sw*w-trace_omega*y)/2;
                        for (int i=0; i<form_count; ++i)
                            if ((su*u+forms[i].trace*x+forms[i].omega_trace*y)%prime == 0)
                                counts[i] += weight;
                    }
            }
        }
        for (const auto c : counts) std::cout << c << ' ';
        std::cout << '\n';
    }
}
