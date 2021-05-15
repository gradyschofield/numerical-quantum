//
// Created by Grady Schofield on 5/13/21.
//

#ifndef ATOM_MRG32K3A_H
#define ATOM_MRG32K3A_H

class Mrg32k3a {
    vector<int64_t> xstate;
    vector<int64_t> ystate;
    static constexpr int64_t m1 = ((int64_t)1 << 31) - 209;
    static constexpr int64_t m2 = ((int64_t)1 << 31) - 22853;
    static constexpr float maxInv = 1.0 / m1;
    static constexpr float maxInvDouble = 1.0 / ((double)m1*m1);

public:

    Mrg32k3a()
            : xstate(3), ystate(3)
    {
        initState();
    }

    void initState() {
        for(auto & y : ystate) {
            y = rand() % m2;
        }
        for(auto & x : xstate) {
            x = rand() % m1;
        }
    }

    void initState(vector<int64_t> const & state, int subsequence) {
        for(int i = 0; i < 3; ++i) {
            xstate[i] = state[subsequence*6+i];
        }
        for(int i = 0; i < 3; ++i) {
            ystate[i] = state[subsequence*6+i+3];
        }
    }

    vector<int64_t> getXState() const {
        return xstate;
    }

    vector<int64_t> getYState() const {
        return ystate;
    }

    static constexpr int64_t getXModulus() {
        return m1;
    }

    static constexpr int64_t getYModulus() {
        return m2;
    }

    uint32_t generate() {
        int64_t x = (xstate[1] * 1403580 - xstate[2] * 810728) % m1;
        if(x < 0) {
            x += m1;
        }
        int64_t y = (ystate[0] * 527612 - ystate[2] * 1370589) % m2;
        if(y < 0) {
            y += m2;
        }
        xstate[2] = xstate[1];
        xstate[1] = xstate[0];
        xstate[0] = x;
        ystate[2] = ystate[1];
        ystate[1] = ystate[0];
        ystate[0] = y;
        int64_t t = (x - y) % m1;
        if(t < 0) {
            t += m1;
        }
        return t;
    }

    float generateFloat() {
        return generate() * maxInv;
    }

    double generateDouble() {
        uint64_t t = generate();
        return (t * m1 + generate()) * maxInvDouble;
    }

    static Matrix3x3 getXMatrix() {
        return Matrix3x3(0, 1, 0, 1403580, 0, 1, -810728, 0, 0);
    }

    static Matrix3x3 getYMatrix() {
        return Matrix3x3(527612, 1, 0, 0, 0, 1, -1370589, 0, 0);
    }

    vector<int64_t> jumpAhead(int64_t totalSamples, int numThreads) {
        vector<int64_t> xstateTmp = getXState();
        vector<int64_t> ystateTmp = getYState();
        Matrix3x3 Ax = getXMatrix();
        Matrix3x3 Ay = getYMatrix();

        int64_t samplesPerThread = totalSamples / numThreads;
        cout << "Samples per thread: " << samplesPerThread << "\n";
        int numSquarings = ceil(log2(samplesPerThread));
        cout << "Skipping ahead " << pow(2, numSquarings) << " by squaring transition matrices " << numSquarings << " times\n";
        timespec t1, t2;
        int N = numSquarings;
        clock_gettime(CLOCK_REALTIME, & t1);
        for(int i = 0; i < N; ++i) {
            Ax.square(getXModulus());
            Ay.square(getYModulus());
        }
        vector<int64_t> state;
        state.reserve(6 * numThreads);
        for(int i = 0; i < numThreads; ++i) {
            for(int j = 0; j < 3; ++j) {
                state.push_back(xstateTmp[j]);
            }
            for(int j = 0; j < 3; ++j) {
                state.push_back(ystateTmp[j]);
            }
            xstateTmp = Ax.multiply(xstateTmp, getXModulus());
            ystateTmp = Ay.multiply(ystateTmp, getYModulus());
        }
        clock_gettime(CLOCK_REALTIME, & t2);
        double time = (t2.tv_sec * 1E9 + t2.tv_nsec - t1.tv_sec * 1E9 - t1.tv_nsec) / 1E9;
        cout << "Time to skip 2^" << N << " for " << numThreads << " threads " << time << "\n";
        return state;
    }
};

#endif //ATOM_MRG32K3A_H
