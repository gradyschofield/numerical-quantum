//
// Created by grady on 5/8/21.
//

#include <iostream>
#include <iomanip>
#include <random>
#include <unistd.h>

#include<cuda.h>
#include<curand_kernel.h>

#include<CommandLine.h>
#include<CudaArray.h>
#include<CudaModule.h>
#include<CudaFunction.h>

using namespace std;

int main(int argc, char ** argv) {
    long N = CommandLine::parseNumPoints(argc, argv, 1E10);

    CudaModule module("FirstIntegral.ptx");
    CudaFunction function = module.getFunction("firstIntegral");

    int numThreads = 64 * 256;

    CudaArray<curandState> curandStates(numThreads);
    CudaArray<double> inTotal(numThreads);
    CudaArray<double> integralTotal(numThreads);
    CudaArray<double> outTotal(numThreads);
    CudaArray<uint64_t> seed(1);
    CudaArray<uint64_t> numSamples(1);

    numSamples[0] = N;
    numSamples.upload();

    random_device r;
    seed[0] = r();
    seed.upload();

    void * args[] = {curandStates.getDevicePtr(),
                     seed.getDevicePtr(),
                     numSamples.getDevicePtr(),
                     inTotal.getDevicePtr(),
                     integralTotal.getDevicePtr(),
                     outTotal.getDevicePtr()};

    cout << "Starting integration" << endl;
    timespec t1, t2;
    clock_gettime(CLOCK_REALTIME, &t1);
    function.run(numThreads, args);
    function.wait();
    clock_gettime(CLOCK_REALTIME, &t2);
    sleep(2);

    inTotal.download();
    integralTotal.download();
    outTotal.download();

    double totalIn = 0, totalIntegral = 0, totalOut = 0;
    for(int i = 0; i < numThreads; ++i) {
        totalIn += inTotal[i];
        totalIntegral += integralTotal[i];
        totalOut += outTotal[i];
    }
    float referenceVolume = pow(4 * M_PI * pow(2,3) / 3, 2);
    cout << "total in: " << totalIn << "\n";
    cout << "total out: " << totalOut << "\n";
    cout << "total integral: " << totalIntegral << "\n";
    float volumeNormalizer = referenceVolume / (totalIn + totalOut);
    double approx = totalIntegral * volumeNormalizer;
    cout << setprecision(16) << "Monte carlo integral: " <<  approx << "\n";
    cout << "Relative error: " <<  fabs(approx-4*M_PI*M_PI)/(4*M_PI*M_PI) << "\n";
    double time = (t2.tv_sec * 1E9 + t2.tv_nsec - t1.tv_sec * 1E9 - t1.tv_nsec)/1E9;
    cout << "time: " << time << "\n";
    cout << "time per sample: " << time / numSamples[0] << "\n";

    return 0;
}
