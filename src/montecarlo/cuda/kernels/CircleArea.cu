#include<curand_kernel.h>

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