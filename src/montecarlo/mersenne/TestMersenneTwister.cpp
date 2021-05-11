//
// Created by Grady Schofield on 5/10/21.
//

#include<iostream>
#include<vector>

using namespace std;

class MersenneTwister {
    /*
     * The naming convention follows the paper http://www.math.sci.hiroshima-u.ac.jp/m-mat/MT/ARTICLES/mt.pdf
     * See the algorithm pseudocode on page 7 and Table 2 page 8 of that paper.
     */
    vector<uint32_t> state;
    int n;
    int m;
    int a;
    int u;
    int s;
    int t;
    int l;
    int b;
    int c;

    int lastPos = 0;
    uint32_t upperMask;
    uint32_t lowerMask;
public:
    MersenneTwister(int n, int m, int r, int a, int u, int s, int t, int l, int b, int c)
        : state(n), n(n), m(m), a(a), u(u), s(s), t(t), l(l), b(b), c(c)
    {
        upperMask = ((uint32_t)~0) << (32-r);
        lowerMask = ((uint32_t)~0) >> r;
    }

    void initState() {
        for(uint32_t & x : state) {
            x = rand();
        }
    }

    uint32_t generate() {
        uint32_t first = state[lastPos];
        uint32_t second = state[lastPos + 1 >= n ? lastPos + 1 - n : lastPos + 1];
        uint32_t mth = state[lastPos + m >= n ? lastPos + m - n : lastPos + m];

        uint32_t t1 = (upperMask & first) | (lowerMask & second);
        uint32_t t2 = mth ^ (t1 >> 1);
        uint32_t x = (t1 & 0x1) == 0 ? t2 : t2 ^ a;
        state[lastPos] = x;
        uint32_t y = x ^ (x >> u);
        y = y ^ ((y << s) & b);
        y = y ^ ((y << t) & c);
        y = y ^ (y >> l);
        lastPos = lastPos + 1 == n ? 0 : lastPos + 1;
        return y;
    }

    float generateFloat() {
        static float maxInv = 1.0f / (float)numeric_limits<uint32_t>::max();
        uint32_t t = generate();
        return t * maxInv;
    }

    double generateDouble() {
        static double maxInv = 1.0 / (double)numeric_limits<uint64_t>::max();
        uint64_t t1 = generate();
        uint64_t t2 = generate();
        return (t1 | t2 << 32) * maxInv;
    }
};

int main(int argc, char ** argv) {
    MersenneTwister mt19937(624, 397, 31, 0x9908b0df, 11, 7, 15, 18, 0x9d2c5680, 0xefc60000);
    mt19937.initState();
    timespec t1, t2;
    uint64_t n = 1E7;
    clock_gettime(CLOCK_REALTIME, &t1);
    for(uint64_t i = 0; i < n; ++i) {
        mt19937.generateFloat();
    }
    clock_gettime(CLOCK_REALTIME, &t2);
    double time = (t2.tv_sec * 1E9 + t2.tv_nsec - t1.tv_sec * 1E9 - t1.tv_nsec)/1E9;
    cout << time / n << " seconds per sample\n";
    cout << 6 * time / n << " seconds per pair of 3 vector\n";
    return 0;
}