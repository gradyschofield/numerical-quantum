//
// Created by grady on 5/6/21.
//

#ifndef CUDA_CUDACONTEXT_H
#define CUDA_CUDACONTEXT_H

#include<atomic>
#include<mutex>
#include<sstream>
#include<iostream>
#include<memory>

#include<cuda.h>

using namespace std;

class CudaContext {
    static mutex contextMutex;
    static shared_ptr<CudaContext> cudaContext;

    CUcontext context = nullptr;

    CudaContext() {
        auto checkError = [](CUresult ret, string function) {
            if(ret != 0) {
                stringstream sstr;
                sstr << "Got CUDA error " << ret << " on " << function << " in CudaContext contructor.";
                throw new runtime_error(sstr.str());
            }
        };

        CUresult ret = cuInit(0);
        checkError(ret, "cuInit");
        CUdevice device;
        ret = cuDeviceGet(&device, 0);
        checkError(ret, "cuDeviceGet");
        ret = cuCtxCreate(&context, 0, device);
        checkError(ret, "cuCtxCreate");
    }

public:

    virtual ~CudaContext() {
        if(context) {
            cout << "Destructor called" << endl;
            cuCtxDestroy(context);
            context = nullptr;
        }
    }

    static shared_ptr<CudaContext> getContext() {
        lock_guard<mutex> lock(contextMutex);
        if(!cudaContext) {
            cudaContext.reset(new CudaContext());
        }
        return cudaContext;
    }
};


#endif //CUDA_CUDACONTEXT_H
