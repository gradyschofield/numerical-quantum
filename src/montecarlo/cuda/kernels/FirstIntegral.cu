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

    double in = 0, out = 0, integral = 0;
    for(unsigned long long i = 0; i < n; ++i) {
        float k_r = 2 * curand_uniform(state);
        float k_theta = M_PI * curand_uniform(state);
        float k_phi = 2 * M_PI * curand_uniform(state);

        float q_r = 2 * curand_uniform(state);
        float q_theta = M_PI * curand_uniform(state);
        float q_phi = 2 * M_PI * curand_uniform(state);

        float volumeProduct = k_r * k_r * sin(k_theta) * q_r * q_r * sin(q_theta);

        float k_x = k_r * cos(k_phi) * sin(k_theta);
        float k_y = k_r * sin(k_phi) * sin(k_theta);
        float k_z = k_r * cos(k_theta);

        float q_x = q_r * cos(q_phi) * sin(q_theta);
        float q_y = q_r * sin(q_phi) * sin(q_theta);
        float q_z = q_r * cos(q_theta);

        float sum_x = k_x + q_x;
        float sum_y = k_y + q_y;
        float sum_z = k_z + q_z;
        float sum_len2 = sum_x * sum_x + sum_y * sum_y + sum_z * sum_z;

        out += (k_r > 1 || sum_len2 > 1) ? volumeProduct : 0;
        integral += (k_r <= 1 && sum_len2 <= 1) ? volumeProduct/(q_r*q_r) : 0;
        in += (k_r <= 1 && sum_len2 <= 1) ? volumeProduct : 0;
    }
    totalInVolume[seq] = in;
    totalOutVolume[seq] = out;
    totalIntegral[seq] = integral;
}