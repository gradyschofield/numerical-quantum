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
        float x = curand_uniform(state) - 0.5f;
        float y = curand_uniform(state) - 0.5f;

        int z = x*x + y*y < 0.25f ? 1 : 0;
        in += z;
        out += 1-z;
    }
    inCount[seq] = in;
    outCount[seq] = out;
}