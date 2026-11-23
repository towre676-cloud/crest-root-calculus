// Exact point counts for the two signed triple factors (0,1,2) at p=11.
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <vector>

using u32 = std::uint32_t;
using i64 = std::int64_t;
constexpr int p = 11;

class Field {
public:
    int degree;
    u32 size = 1, order;
    std::vector<u32> exponential, logarithm;
    std::array<std::vector<std::uint8_t>,2> cubics;
    explicit Field(int n) : degree(n) {
        if (n<1 || n>8) throw std::runtime_error("Extension must be between one and eight");
        for (int i=0;i<n;++i) size*=p;
        order=size-1;
        // Coefficients at X^0, X^1, X^2 of the monic primitive modulus.
        static const int moduli[9][3] = {{0,0,0},{0,0,0},{2,4,0},{3,5,0},
                                              {2,1,0},{9,0,2},{2,2,4},{3,3,0},{2,1,4}};
        exponential.resize(size);
        logarithm.assign(size,order);
        u32 top=size/p;
        constexpr u32 radix=p*p*p;
        std::array<std::array<u32,radix>,p> reduced{};
        for(int lead=0;lead<p;++lead) for(u32 lo=0;lo<radix;++lo) {
            u32 v=lo,encoded=0,power=1;
            for(int j=0;j<3;++j) {
                encoded+=power*((int(v%p)+p*p-lead*moduli[n][j])%p);
                power*=p; v/=p;
            }
            reduced[lead][lo]=encoded;
        }
        u32 x=1;
        for(u32 e=0;e<order;++e) {
            if(logarithm[x]!=order) throw std::runtime_error("Primitive orbit repeated");
            exponential[e]=x; logarithm[x]=e;
            if(n==1) x=x*2%p;
            else {
                u32 lead=x/top, shifted=(x%top)*p;
                x=shifted-shifted%radix+reduced[lead][shifted%radix];
            }
        }
        if(x!=1 || logarithm[0]!=order) throw std::runtime_error("Invalid field tables");
        exponential[order]=0;
        for(auto& t:cubics) t.assign(size,0);
        for(auto& t:cubics) t[order]=1; // x=0 in both cubic histograms
        for(u32 e=0;e<order;++e) {
            const u32 sq=exponential[(2*e)%order];
            for(int c=0;c<2;++c) {
                // x^3+c*x = x*(x^2+c), with c either 1 or the generator.
                u32 y=logarithm[add_class(sq,c)];
                u32 image=y==order?order:(e+y)%order;
                ++cubics[c][image];
            }
        }
    }
    u32 translate(u32 x,int c) const { return x-x%p+(x%p+c)%p; }
    u32 add_class(u32 x,int c) const {
        if(c==0 || degree==1) return translate(x,c==0?1:2);
        u32 digit=(x/p)%p;
        return x-digit*p+((digit+1)%p)*p;
    }
    u32 residue(i64 a) const {
        i64 v=a%order;
        return static_cast<u32>(v<0?v+order:v);
    }
    u32 shift(u32 a,i64 b) const { return a==order?order:residue(i64(a)+b); }
    u32 multiply(u32 a,u32 b) const {
        return a==order || b==order ? order : (a+b)%order;
    }
    int roots(u32 alpha,u32 beta) const {
        if(alpha==order) {
            if(beta==order) return 1;
            int g=std::gcd(3u,order);
            return (beta%g==0)?g:0;
        }
        int odd=alpha%2;
        return cubics[odd][shift(beta,-3LL*((alpha-odd)/2))];
    }
};

struct Counts {
    Field field;
    std::array<i64,2> full{}, signed_sum{};
    i64 seen=0;
    std::vector<int> digits;
    explicit Counts(int d):field(d),digits(d+1){}
    void point(u32 e,int orbit) {
        const auto& F=field;
        u32 z=F.exponential[e];
        u32 u=F.logarithm[p-4];
        for(int j=0;j<3;++j)
            u=F.logarithm[F.translate(F.exponential[F.multiply(u,e)],1)];
        seen+=orbit;
        for(int twist=0;twist<2;++twist) {
            if(u!=F.order && u%2!=static_cast<u32>(twist)) continue;
            i64 weight=i64(orbit)*(u==F.order?1:2);
            u32 beta=u==F.order?F.order:F.residue(i64((u-twist)/2)-twist);
            int all=1, reduced=1;
            for(int i=0;i<3;++i) {
                u32 a=F.shift(F.logarithm[F.translate(z,i)],i64(F.logarithm[3])-twist);
                int v=F.roots(a,beta);
                all*=v; reduced*=v-1;
            }
            full[twist]+=weight*all;
            signed_sum[twist]+=weight*reduced;
        }
    }
    void necklaces(int t,int period) {
        if(t>field.degree) {
            if(field.degree%period==0) {
                u32 e=0;
                for(int i=1;i<=field.degree;++i) e=p*e+digits[i];
                if(e<field.order) point(e,period);
            }
            return;
        }
        digits[t]=digits[t-period];
        necklaces(t+1,period);
        for(int j=digits[t-period]+1;j<p;++j) {
            digits[t]=j; necklaces(t+1,t);
        }
    }
    void compute() {
        point(field.order,1); // z=0
        necklaces(1,1);
        if(seen!=field.size) throw std::runtime_error("Frobenius orbit sizes failed");
        for(int twist=0;twist<2;++twist) {
            int all=1, reduced=1;
            for(int i=0;i<3;++i) {
                const int a=(11-(12*i+1)%11)%11; // 12 = 1 modulo 11
                int chi=field.logarithm[a]%2==0 ? 1 : -1;
                int sheets=2+(twist==0?chi:-chi);
                all*=sheets; reduced*=sheets-1;
            }
            full[twist]+=all; signed_sum[twist]+=reduced;
        }
        std::cout<<field.degree<<' '<<field.size<<' '<<full[0]<<' '<<full[1]<<' '
                 <<-(signed_sum[0]+signed_sum[1])/2<<' '
                 <<-(signed_sum[0]-signed_sum[1])/2<<'\n';
    }
};

int main(int argc,char** argv) {
    if(argc!=2) {std::cerr<<"Supply an extension degree from 1 through 8\n";return 2;}
    Counts(std::stoi(argv[1])).compute();
}
