//
// Created by Grady Schofield on 5/13/21.
//

#include<iostream>
#include<vector>
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

    int64_t numSamples = 10E12;
    int64_t numThreads = 1024;

    vector<int64_t> threadStates = mrg32K3A.jumpAhead(numSamples, numThreads);

    return 0;
}