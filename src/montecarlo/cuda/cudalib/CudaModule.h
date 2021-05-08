//
// Created by grady on 5/6/21.
//

#ifndef CUDA_CUDAMODULE_H
#define CUDA_CUDAMODULE_H

#include<cuda.h>

#include<CudaContext.h>
#include<CudaFunction.h>

class CudaModule {
    CUmodule module;
    shared_ptr<CudaContext> cudaContext;

public:

    CudaModule(string filename) {
        cudaContext = CudaContext::getContext();
        CUresult ret = cuModuleLoad(&module, filename.c_str());
        if(ret != 0) {
            stringstream sstr;
            sstr << "Got CUDA error " << ret << " on cuModuleLoad in CudaModule contructor.";
            throw runtime_error(sstr.str());
        }
    }

    CudaFunction getFunction(string name) {
        CUfunction function;
        CUresult ret = cuModuleGetFunction(&function, module, name.c_str());
        if(ret != 0) {
            stringstream sstr;
            sstr << "Got CUDA error " << ret << " on cuModuleGetFunction in CudaModule::getFunction.";
            throw runtime_error(sstr.str());
        }
        return CudaFunction(function);
    }

    virtual ~CudaModule() {
        if(module) {
            cuModuleUnload(module);
            module = nullptr;
            cout << "Module unloaded" << endl;
        }
    }
};

#endif //CUDA_CUDAMODULE_H
