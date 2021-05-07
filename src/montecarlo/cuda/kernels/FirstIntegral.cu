#include<curand_kernel.h>

#include<BigSum.cuh>

extern "C" __global__ void firstIntegral(curandState * states,
                                         unsigned long long * seed,
                                         unsigned long long * numSamples,
                                         double * totalInVolume,
                                         double * totalIntegral,
                                         double * totalOutVolume) {

    int numThreads = blockDim.x * gridDim.x;
    unsigned long long n = *numSamples / numThreads;
    unsigned long long seq = threadIdx.x + blockIdx.x * blockDim.x;

    curandState * state = &states[seq];
    curand_init(*seed, seq, 0, state);

    BigSum in, out, integral;
    for(unsigned long long i = 0; i < n; ++i) {

    }
    totalInVolume[seq] = in;
    totalOutVolume[seq] = out;
    totalIntegral[seq] = integral;
}