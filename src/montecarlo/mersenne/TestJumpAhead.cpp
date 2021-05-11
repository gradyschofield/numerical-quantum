//
// Created by Grady Schofield on 5/10/21.
//

#include<iostream>
#include<unordered_set>
#include<vector>
#include<tuple>

using namespace std;

/*
 * This program implements the fast jump ahead idea in
 * http://www.iro.umontreal.ca/~lecuyer/myftp/papers/jumpf2.pdf .
 * It also implements the slower straightforward way to jump ahead as a reference
 * to compare the fast code for correctness.
 */

class MatrixElement {
    int row;
    int column;

public:

    MatrixElement() = delete;

    MatrixElement(MatrixElement && e)
        : row(e.getRow()), column(e.getColumn())
    {
    }

    MatrixElement(MatrixElement const & e)
        : row(e.getRow()), column(e.getColumn())
    {
    }

    MatrixElement(int row, int column)
        : row(row), column(column)
    {
    }

    int getRow() const {
        return row;
    }

    int getColumn() const {
        return column;
    }

    bool operator==(MatrixElement const & e) const {
        return row == e.getRow() && column == e.getColumn();
    }
};

namespace std {
    template<> struct hash<MatrixElement>
    {
        size_t operator()(MatrixElement const & elem) const noexcept {
            return elem.getRow() ^ (elem.getColumn() << 1);
        }
    };
}

/*
 * The matrix is stored in column major ordering.  The first word is the first 32 rows of column 1.
 */
class DenseF2Matrix {
    vector<uint32_t> elements;
    int dim;
    int wordsPerColumn;
public:

    DenseF2Matrix(int dim)
        : dim(dim), wordsPerColumn(dim/32), elements(wordsPerColumn*dim)
    {
    }

    DenseF2Matrix(vector<uint32_t> && elements) {
        if((int)sqrt(elements.size()) != sqrt(elements.size())) {
            throw runtime_error("Tried to construct a nonsquare DenseF2Matrix");
        }
        dim = sqrt(elements.size());
        DenseF2Matrix::elements = move(elements);
    }

    DenseF2Matrix square() const {
    }

    void setElement(int row, int column) {
        int columnWord = column / 32;
        int columnBit = column % 32;
        elements[row * wordsPerColumn + columnWord] |= 1 << columnBit;
    }

    uint32_t getElement(int row, int column) const {
        int columnWord = column / 32;
        int columnBit = column % 32;
        return (elements[row * wordsPerColumn + columnWord] & (1 << columnBit)) == 0 ? 0 : 1;
    }
};

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
        int maxColumn = 0;
        for(auto & p : elements) {
            maxRow = max(maxRow, p.getRow());
            maxColumn = max(maxColumn, p.getColumn());
        }

        if(maxRow != maxColumn) {
            throw runtime_error("Max row didn't match max column in SparseMatrix::bake");
        }

        if(maxRow % w != 0) {
            throw runtime_error("Max row is not a multiple of w");
        }

        if(maxColumn % w != 0) {
            throw runtime_error("Max column is not a multiple of w");
        }

        DenseF2Matrix matrix(maxRow);
        for(auto & p : elements) {
            matrix.setElement(p.getRow(), p.getColumn());
        }
        return matrix;
    }
};

int main(int argc, char ** argv) {
    int n = 624;
    int w = 32;
    int r = 31;
    uint32_t a = 0x9908b0df;
    int rowOffset = 0;
    int columnOffset = 0;
    SparseMatrix A;
    for(int z = 0; z < n; ++z) {
        /*
         * x_k = x_k+m + x_1^u | x_2^l)
         * 0 0 0 ... 0 0 a31
         * 1 0 0 ... 0 0 a30
         * 0 1 0 ... 0 0 a29
         * 0 0 1 ... 0 0 a28
         * . . . .   . .  .
         * . . .  .  . .  .
         * . . .   . . .  .
         * 0 0 0 ... 1 0  a1
         * 0 0 0 ... 0 1  a0
         */

        /*
         * Handle the upper w-r bits of x_k
         */
        for(int row = rowOffset + 1; row < rowOffset + w - r; ++row) {
            int col = row - 1;
            A.insert(row, col);
        }

        /*
         * Handle the lower r bits of x_(k+1)
         */
        for(int row = rowOffset + w - r; row < rowOffset + w; ++row) {
            int col = row + w - 1;
            A.insert(row, col);
        }

        for(int i = 0; i < w; ++i) {
            int row = rowOffset + i;
            int col = columnOffset + 2*w;
            // The low order bit of 'a' is in the lower left corner of the matrix, hence a shift of w-i
            if(a & (1 << (w-i))) {
                A.insert(row, col);
            }
        }

        columnOffset += w;
        rowOffset += w;
    }

    DenseF2Matrix denseMatrix = A.bake(w);

    cout << "Num nonzero elements in A " << A.getNumElements() << "\n";
    return 0;
}