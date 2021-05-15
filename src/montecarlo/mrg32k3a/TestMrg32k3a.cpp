//
// Created by Grady Schofield on 5/13/21.
//

#include<iostream>
#include<vector>
#include<thread>
#include<cmath>

#include<Matrix3x3.h>
#include<Mrg32k3a.h>

using namespace std;

int main(int argc, char ** argv) {
    Mrg32k3a mrg32K3A;
    vector<int64_t> xstateTmp = mrg32K3A.getXState();
    vector<int64_t> ystateTmp = mrg32K3A.getYState();
    Matrix3x3 Ax = Mrg32k3a::getXMatrix();
    Matrix3x3 Ay = Mrg32k3a::getYMatrix();

    for(int i = 0; i < 10; ++i) {
        xstateTmp = Ax.multiply(xstateTmp, Mrg32k3a::getXModulus());
        ystateTmp = Ay.multiply(ystateTmp, Mrg32k3a::getYModulus());
        cout << "generate: " << mrg32K3A.generate() << " ";
        int64_t t = (xstateTmp[0] - ystateTmp[0]) % Mrg32k3a::getXModulus();
        if(t < 0) {
            t += Mrg32k3a::getXModulus();
        }
        cout << "manual: " << t << "\n";
    }

    cout << "\n";

    int64_t numSamples = 10E12;
    int64_t numThreads = 1024;

    vector<int64_t> threadStates = mrg32K3A.jumpAhead(numSamples, numThreads);

    cout << "\n";

    numSamples = 1E10;
    numThreads = thread::hardware_concurrency();
    threadStates = mrg32K3A.jumpAhead(numSamples, numThreads);

    cout << "\n";

    timespec t1, t2;
    double t = 0;
    long N = numSamples / numThreads;
    clock_gettime(CLOCK_REALTIME, &t1);
    vector<double> accum(numThreads);
    vector<thread> threads;
    cout << "Computing " << (double)numSamples << " numbers with " << numThreads << " threads...\n";
    for(int threadIdx = 0; threadIdx < numThreads; ++threadIdx) {
        threads.emplace_back([threadIdx,N,&threadStates,&accum](){
            Mrg32k3a generator;
            generator.initState(threadStates, threadIdx);
            double t = 0;
            for(long i = 0; i < N; ++i) {
                t += generator.generateFloat();
            }
            accum[threadIdx] = t;
        });
    }
    for(int i = 0; i < numThreads; ++i) {
        threads[i].join();
        t += accum[i];
    }
    clock_gettime(CLOCK_REALTIME, &t2);
    cout << t << "\n";
    double time = (t2.tv_sec * 1E9 + t2.tv_nsec - t1.tv_sec * 1E9 - t1.tv_nsec)/1E9;
    cout << time / (N * numThreads) << " seconds per sample\n";

    return 0;
}