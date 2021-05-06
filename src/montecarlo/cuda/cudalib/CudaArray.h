//
// Created by grady on 5/6/21.
//

#ifndef CUDA_CUDAARRAY_H
#define CUDA_CUDAARRAY_H

#include<vector>
#include<memory>
#include<cuda.h>

using namespace std;

template<typename T>
class CudaArray {
    vector<T> hostArray;
    unique_ptr<CUdeviceptr> devicePtr;

public:

    CudaArray(int size)
        : hostArray(size), devicePtr(make_unique<CUdeviceptr>())
    {
        cuMemAlloc(devicePtr.get(), sizeof(T) * size);
        upload();
    }

    CUdeviceptr * getDevicePtr() {
        return devicePtr.get();
    }

    T & operator[](int i) {
        return hostArray[i];
    }

    const T & operator[](int i) const {
        return hostArray[i];
    }

    int upload() const {
        int ret = cuMemcpyHtoD(*devicePtr, hostArray.data(), sizeof(T) * hostArray.size());
        if(ret != 0) {
            cout << "cuMemcpyHtoD error " << ret << "\n";
        }
        return ret;
    }

    int download() {
        int ret = cuMemcpyDtoH(hostArray.data(), *devicePtr, sizeof(T) * hostArray.size());
        if(ret != 0) {
            cout << "cuMemcpyDtoH error " << ret << "\n";
        }
        return ret;
    }

    void free() {
        if(devicePtr) {
            cuMemFree(*devicePtr);
        }
        hostArray = vector<T>();
        devicePtr.reset(nullptr);
    }

    virtual ~CudaArray() {
        if(devicePtr) {
            cuMemFree(*devicePtr);
        }
    }
};

#endif //CUDA_CUDAARRAY_H
