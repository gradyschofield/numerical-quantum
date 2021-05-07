#include<curand_kernel.h>

class BigSum {
    float positiveAccumulators[25];
    float negativeAccumulators[25];
    float orderLimits[25];
    int binOffset;
    int orders;

public:

    BigSum() {
        float minBin = 1E-5;
        orders = 25;
        binOffset = round(log10(minBin));
        double lowestLimit = pow(10, binOffset);
        for (int i = 0; i < orders; ++i) {
            orderLimits[i] = lowestLimit * pow(10, i);
            positiveAccumulators[i] = 0;
            negativeAccumulators[i] = 0;
        }
    }
    void add(float x) {
        int bin = min(orders-1, (int)max(0.0, round(log10(fabs(x))) - binOffset));
        float positiveInc = max(0.0f, x);
        float negativeInc = min(0.0f, x);
        positiveAccumulators[bin] += positiveInc;
        negativeAccumulators[bin] += negativeInc;
    }

    operator double() const {
        double positiveSum = 0;
        double negativeSum = 0;
        for(int i = 0; i < orders; ++i) {
            positiveSum += positiveAccumulators[i];
        }
        for(int i = 0; i < orders; ++i) {
            negativeSum += negativeAccumulators[i];
        }
        return positiveSum + negativeSum;
    }

    void operator+=(double x) {
        add(x);
    }

    void operator-=(double x) {
        add(-x);
    }

    void carryAccumulators() {
        for (int bin = 0; bin < orders - 1; ++bin) {
            float setValue = positiveAccumulators[bin] > orderLimits[bin] ? 0 : positiveAccumulators[bin];
            float carryValue = positiveAccumulators[bin] > orderLimits[bin] ? positiveAccumulators[bin] : 0;
            positiveAccumulators[bin + 1] += carryValue;
            positiveAccumulators[bin] = setValue;
        }
        for (int bin = 0; bin < orders - 1; ++bin) {
            float setValue = negativeAccumulators[bin] < -orderLimits[bin] ? 0 : negativeAccumulators[bin];
            float carryValue = negativeAccumulators[bin] < -orderLimits[bin] ? negativeAccumulators[bin] : 0;
            negativeAccumulators[bin + 1] += carryValue;
            negativeAccumulators[bin] = setValue;
        }
    }

};

extern "C" __global__ void integrate(curandState * states,
                                     unsigned long long * seed,
                                     unsigned long long * numSamples,
                                     unsigned long long * inCount,
                                     unsigned long long * outCount) {
    int numThreads = blockDim.x * gridDim.x;
    unsigned long long n = *numSamples / numThreads;
    unsigned long long seq = threadIdx.x + blockIdx.x * blockDim.x;

    curandState * state = &states[seq];
    curand_init(*seed, seq, 0, state);

    unsigned long long in = 0;
    unsigned long long out = 0;
    for(unsigned long long i = 0; i < n; ++i) {
        float x = 5 * curand_uniform(state) - 2;
        float y = 5 * curand_uniform(state) - 2;

        float d1 = x*x + y*y;
        float d2 = (x-1)*(x-1) + y*y;
        float d3 = (x-1)*(x-1) + (y-1)*(y-1);
        float d4 = x*x + (y-1)*(y-1);

        int z = d1 < 4 && d2 < 4 && d3 < 4 && d4 < 4 ? 1 : 0;
        in += z;
        out += 1-z;
    }
    inCount[seq] = in;
    outCount[seq] = out;
}