//
// Created by Grady Schofield on 4/25/21.
//

#ifndef ATOM_MATHUTIL_H
#define ATOM_MATHUTIL_H

#include<time.h>

template<typename T>
T square(T x) {
    return x * x;
}

template<typename T>
T cube(T x) {
    return x * x * x;
}

class Time {
public:
    static timespec startTimer() {
        timespec t;
        clock_gettime(CLOCK_REALTIME, &t);
        return t;
    }

    static double stopTimer(timespec const &start) {
        timespec t;
        clock_gettime(CLOCK_REALTIME, &t);
        return ((t.tv_sec - start.tv_sec) * 1E9 + t.tv_nsec - start.tv_nsec)/1E9;
    }
};

#endif //ATOM_MATHUTIL_H
