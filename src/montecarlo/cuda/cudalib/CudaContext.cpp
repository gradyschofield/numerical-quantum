//
// Created by grady on 5/6/21.
//

#include<atomic>
#include<memory>

#include<cuda.h>

#include "CudaContext.h"

using namespace std;

shared_ptr<CudaContext> CudaContext::cudaContext;
mutex CudaContext::contextMutex;
