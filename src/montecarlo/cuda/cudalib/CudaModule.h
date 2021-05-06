//
// Created by grady on 5/6/21.
//

#ifndef CUDA_CUDAMODULE_H
#define CUDA_CUDAMODULE_H

#include<cuda.h>

class CudaModule {
    CUModule module;

public:
    virtual ~CudaModule() {
        cuModuleUnload(module);
    }
};

#endif //CUDA_CUDAMODULE_H
