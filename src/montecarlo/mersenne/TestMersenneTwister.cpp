//
// Created by Grady Schofield on 5/10/21.
//

#include<iostream>
#include<vector>

#include<MersenneTwister.h>

using namespace std;

int main(int argc, char ** argv) {
    MersenneTwister mt19937 = MersenneTwister::createMT19937();
    timespec t1, t2;
    uint64_t n = 1E7;
    clock_gettime(CLOCK_REALTIME, &t1);
    for(uint64_t i = 0; i < n; ++i) {
        mt19937.generateFloat();
    }
    clock_gettime(CLOCK_REALTIME, &t2);
    double time = (t2.tv_sec * 1E9 + t2.tv_nsec - t1.tv_sec * 1E9 - t1.tv_nsec)/1E9;
    cout << time / n << " seconds per sample\n";
    cout << 6 * time / n << " seconds per pair of 3 vector\n";
    return 0;
}