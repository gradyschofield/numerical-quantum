//
// Created by grady on 5/6/21.
//

#ifndef CUDA_CUDAFUNCTION_H
#define CUDA_CUDAFUNCTION_H

#include<cuda.h>

class CudaFunction {
    CUfunction function;
public:
    CudaFunction(CUfunction function)
        : function(function) {
    }

    void run(int numThreads, void * args[]) {
        int numXBlocks = numThreads / 64;
        CUresult ret = cuLaunchKernel(function,
                                      numXBlocks, 1, 1,
                                      64, 1, 1,
                                      0, 0, args, 0);
        if(ret != 0) {
            stringstream sstr;
            sstr << "Got CUDA error " << ret << " on cuLaunchKernel in CudaFunction::run.";
            throw runtime_error(sstr.str());
        }
    }

    void run(int numXBlocks, int numYBlocks, int numZBlocks,
             int numXThreads, int numYThreads, int numZThreads,
             CUstream stream, void * args[]) {
        CUresult ret = cuLaunchKernel(function,
                                      numXBlocks, numYBlocks, numZBlocks,
                                      numXThreads, numYThreads, numZThreads,
                                      0, stream, args, 0);
        if(ret != 0) {
            stringstream sstr;
            sstr << "Got CUDA error " << ret << " on cuLaunchKernel in CudaFunction::run.";
            throw runtime_error(sstr.str());
        }
    }

    void wait() {
        cuStreamSynchronize(0);
    }
};

#endif //CUDA_CUDAFUNCTION_H
