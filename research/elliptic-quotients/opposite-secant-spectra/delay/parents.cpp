#include <algorithm>
#include <cstdlib>
#include <iostream>

using i64 = long long;

// Q = a*x*x+b*y*y+c*z*z+2*(u*x*y+v*x*z+w*y*z).
// The ternary Minkowski inequalities give abc <= 2*det(Q).
// The signs of the basis vectors can always make u and v nonnegative.
int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "usage: parents MAX_CHILD_DETERMINANT\n";
        return 2;
    }
    const i64 limit = std::atoll(argv[1]) / 9;
    if (limit < 1 || limit > 1000000) {
        std::cerr << "parent determinant limit must be between 1 and 1000000\n";
        return 2;
    }
    for (i64 a = 1; a*a*a <= 2*limit; ++a) {
        if (a%4 > 1) continue;
        for (i64 b = a; a*b*b <= 2*limit; ++b) {
            if (b%4 > 1) continue;
            for (i64 u = (a*b)%2; 2*u <= a; u += 2) {
                const i64 minor = a*b-u*u;
                for (i64 v = 0; 2*v <= a; ++v) {
                    for (i64 w = -b/2; 2*w <= b; ++w) {
                        if ((u == 0 || v == 0) && w < 0) continue;
                        // Test all sign choices of x,y,z in {1,-1}.
                        if (2*std::min({u+v+w,u-v-w,-u+v-w,-u-v+w}) < -a-b) continue;
                        const i64 shift = 2*u*v*w-a*w*w-b*v*v;
                        const i64 largest_c = std::min(2*limit/(a*b), (limit-shift)/minor);
                        for (i64 c = b; c <= largest_c; ++c) {
                            if (c%4 > 1 || (v-a*c)%2 != 0 || (w-b*c)%2 != 0) continue;
                            const i64 det = minor*c+shift;
                            // Norm one is internal on the conic exactly for determinant 1 mod 3.
                            if (det <= 0 || det%3 != 1 || a*b*c > 2*det) continue;
                            std::cout << det << ' ' << a << ' ' << b << ' ' << c << ' '
                                      << u << ' ' << v << ' ' << w << '\n';
                        }
                    }
                }
            }
        }
    }
}
