//
// Created by Grady Schofield on 5/12/21.
//

#ifndef ATOM_SPARSEMATRIX_H
#define ATOM_SPARSEMATRIX_H

#include<iostream>
#include<unordered_set>

using namespace std;

#include<DenseF2Matrix.h>
#include<MatrixElement.h>

/*
 * This class will help us build up the recurrence matrix for a Mersenne twister generator
 * using a code that is fairly easy to read.
 */
class SparseMatrix {
    unordered_set<MatrixElement> elements;

public:
    void insert(int row, int col) {
        elements.emplace(row, col);
    }

    int getNumElements() const {
        return elements.size();
    }

    DenseF2Matrix bake(int w) const {
        int maxRow = 0;
        for(auto & p : elements) {
            maxRow = max(maxRow, p.getRow());
        }
        int dim = maxRow + 1;

        if(dim % w != 0) {
            stringstream sstr;
            sstr << "Dimension is not a multiple of w. Dimension : " << dim << " w: " << w;
            throw runtime_error(sstr.str());
        }

        cout << "Max row " << dim << "\n";
        DenseF2Matrix matrix(dim);
        for(auto & p : elements) {
            matrix.setElement(p.getRow(), p.getColumn());
        }
        return matrix;
    }
};

#endif //ATOM_SPARSEMATRIX_H
