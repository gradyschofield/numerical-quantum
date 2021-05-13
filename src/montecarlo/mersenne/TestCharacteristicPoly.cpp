//
// Created by Grady Schofield on 5/13/21.
//

#include<iostream>

using namespace std;

#include<DenseF2Matrix.h>

void printBits(vector<uint32_t> const & v) {
    for(uint32_t x : v) {
        for(int i = 31; i >= 0; --i) {
            cout << ((x & 1 << i) == 0 ? 0 : 1);
        }
        cout << "\n";
    }
}

int main(int argc, char ** argv) {
    vector<uint32_t> poly;
    poly.push_back(0x80000000);
    poly.push_back(0);
    printBits(poly);
    DenseF2Matrix::polyProduct(true, poly);
    cout << "after:\n";
    printBits(poly);

    return 0;
}