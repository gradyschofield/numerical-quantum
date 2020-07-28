//
// Created by grady on 6/22/20.
//

#include<iostream>
#include "Complex.h"

using namespace std;

int main(int argc, char ** argv) {
    Complex<double> a(2.0, -1.0), b(3.0, 1.0);
    cout << a / b << "\n";
    auto c = a/b;
    cout << c * b << "\n";
    cout << exp(c) << "\n";
    return 0;
}
