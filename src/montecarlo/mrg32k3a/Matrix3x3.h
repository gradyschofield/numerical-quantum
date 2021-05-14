//
// Created by Grady Schofield on 5/13/21.
//

#ifndef ATOM_MATRIX3X3_H
#define ATOM_MATRIX3X3_H

#include<vector>

using namespace std;

class Matrix3x3 {
    vector<int64_t> elements;

public:

    Matrix3x3()
            : elements(9)
    {
    }

    Matrix3x3(
            int64_t m11, int64_t m21, int64_t m31,
            int64_t m12, int64_t m22, int64_t m32,
            int64_t m13, int64_t m23, int64_t m33)
            : elements({m11, m21, m31, m12, m22, m32, m13, m23, m33})
    {
    }

    void square(int64_t modulus) {
        Matrix3x3 ret;
        for(int row = 0; row < 3; ++row) {
            for(int column = 0; column < 3; ++column) {
                int64_t dot = 0;
                for(int k = 0; k < 3; ++k) {
                    int64_t term = getElement(row, k) * getElement(k, column) % modulus;
                    if(term < 0) {
                        term += modulus;
                    }
                    dot += term;
                    dot %= modulus;
                    if(dot < 0) {
                        dot += modulus;
                    }
                }
                ret.setElement(row, column, dot);
            }
        }
        *this = ret;
    }

    vector<int64_t> multiply(vector<int64_t> const & v, int64_t modulus) {
        vector<int64_t> ret(3);
        for(int row = 0; row < 3; ++row) {
            for(int column = 0; column < 3; ++column) {
                ret[row] += getElement(row, column) * v[column];
            }
            ret[row] %= modulus;
            if(ret[row] < 0) {
                ret[row] += modulus;
            }
        }
        return ret;
    }

    int64_t getElement(int row, int column) const {
        return elements[column * 3 + row];
    }

    void setElement(int row, int column, int64_t e) {
        elements[column * 3 + row] = e;
    }
};

#endif //ATOM_MATRIX3X3_H
