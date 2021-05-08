#include <iostream>
#include <iomanip>

#include<cuda.h>
#include<curand_kernel.h>

#include<CudaArray.h>
#include<CudaModule.h>
#include<CudaFunction.h>

using namespace std;

int main() {

    CudaModule module("CircleArea.ptx");
    CudaFunction function = module.getFunction("integrate");

    int numThreads = 64 * 256;

    CudaArray<curandState> curandStates(numThreads);
    CudaArray<uint64_t> inCount(numThreads);
    CudaArray<uint64_t> outCount(numThreads);
    CudaArray<uint64_t> seed(1);
    CudaArray<uint64_t> numSamples(1);

    numSamples[0] = 1E12;
    numSamples.upload();

    seed[0] = 846123423;
    seed.upload();

    void * args[] = {curandStates.getDevicePtr(),
                     seed.getDevicePtr(),
                     numSamples.getDevicePtr(),
                     inCount.getDevicePtr(),
                     outCount.getDevicePtr()};

    timespec t1, t2;
    clock_gettime(CLOCK_REALTIME, &t1);
    function.run(numThreads, args);
    function.wait();
    clock_gettime(CLOCK_REALTIME, &t2);

    inCount.download();
    outCount.download();

    uint64_t totalIn = 0, totalOut = 0;
    for(int i = 0; i < numThreads; ++i) {
        totalIn += inCount[i];
        totalOut += outCount[i];
    }
    cout << "total in: " << totalIn << "\n";
    cout << "total out: " << totalOut << "\n";
    double approx = 25.0 * totalIn / (double)(totalIn + totalOut);
    cout << setprecision(16) << approx << "\n";
    cout << "cpu approx: " << 12.5664 << "\n";
    double time = (t2.tv_sec * 1E9 + t2.tv_nsec - t1.tv_sec * 1E9 - t1.tv_nsec)/1E9;
    cout << "time: " << time << "\n";
    cout << "time per sample: " << time / numSamples[0] << "\n";

    return 0;
}
