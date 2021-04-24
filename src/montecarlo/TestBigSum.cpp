//
// Created by grady on 4/24/21.
//

#include<iostream>
#include<iomanip>
#include<time.h>

#include<BigSum.h>

using namespace std;

int main(int argc, char ** argv) {
    BigSum sum1(16, 1E-5);
    double sum = -1;
    sum1.add(-1);
    timespec t1, t2;
    clock_gettime(CLOCK_REALTIME, &t1);
    for(long i = 0; i < 4E8; ++i) {
        sum1 += -1E-17;
        sum += -1E-17;
    }
    clock_gettime(CLOCK_REALTIME, &t2);
    cout << "time: " << ((t2.tv_sec - t1.tv_sec) * 1E9 + t2.tv_nsec - t1.tv_nsec)/1E9 << " (sec)\n";
    sum1.print();
    cout << "double sum: " << setprecision(16) << sum << setprecision(5) << "\n";
    return 0;
}