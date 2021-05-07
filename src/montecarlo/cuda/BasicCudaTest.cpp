#include <iostream>
#include <iomanip>
#include <vector>

#include<cuda.h>
#include<curand_kernel.h>

#include<CudaArray.h>

using namespace std;

int main() {
    CUresult ret;
    ret = cuInit(0);
    cout << "cu init: " << ret << "\n";

    int devCount;
    cuDeviceGetCount(&devCount);
    cout << "device count: " << devCount << "\n";

    CUdevice device;
    ret = cuDeviceGet(&device, 0);
    cout << "device get: " << ret << "\n";

    CUcontext context;
    ret = cuCtxCreate(&context, 0, device);
    cout << "ctx create: " << ret << "\n";

    CUmodule module;
    ret = cuModuleLoad(&module, "CircleArea.ptx");
    cout << "module load: " << ret << "\n";

    CUfunction function;
    ret = cuModuleGetFunction(&function, module, "integrate");
    cout << "get function: " << ret << "\n";

    /*
    CudaModule module("random.ptx");
    CudaFunction function = module.getFunction("initRandom");
     */

    int numThreadsPerBlock = 64;
    int numBlocks = 256;
    int numThreads = numBlocks * numThreadsPerBlock;

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
    /*
    function.run(numBlocks, 1, 1,
                 numThreadsPerBlock, 1, 1,
                 0, args);
    function.wait();
     */
    ret = cuLaunchKernel(function,
                   numBlocks, 1, 1,
                   numThreadsPerBlock, 1, 1,
                   0, 0, args, 0);
    cuStreamSynchronize(0);
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

    ret = cuModuleUnload(module);
    cout << "module unload: " << ret << "\n";

    ret = cuCtxDestroy(context);
    cout << "ctx destroy: " << ret << "\n";

    return 0;
}
