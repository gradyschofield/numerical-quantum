//
// Created by Grady Schofield on 5/10/21.
//

#include<iostream>
#include<sstream>
#include<unordered_set>
#include<vector>
#include<thread>
#include<tuple>

#include<MersenneTwister.h>

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
    int dim;
    int wordsPerColumn;
    vector<uint32_t> elements;
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
        DenseF2Matrix ret(dim);
        vector<thread> threads;
        int numThreads = thread::hardware_concurrency();
        int stop = 0;
        for(int i = 0; i < numThreads; ++i) {
            int start = stop;
            stop = start + (dim / numThreads) + (i < dim % numThreads ? 1 : 0);
            threads.emplace_back([start, stop, &ret, this]() {
                for (int column = start; column < stop; ++column) {
                    ret.setColumn(column, multiply(getColumn(column)));
                }
            });
        }
        for(thread & t : threads) {
            t.join();
        }
        return ret;
    }

    vector<uint32_t> getColumn(int column) const {
        vector<uint32_t> ret(wordsPerColumn);
        for(int row = 0; row < wordsPerColumn; ++row) {
            ret[row] = elements[column * wordsPerColumn + row];
        }
        return ret;
    }

    void setColumn(int column, vector<uint32_t> const & v) {
        for(int row = 0; row < wordsPerColumn; ++row) {
            elements[column * wordsPerColumn + row] = v[row];
        }
    }

    vector<uint32_t> multiply(vector<uint32_t> const & v) const {
        vector<uint32_t> ret(wordsPerColumn);
        for(int column = 0; column < dim; ++column) {
            int columnWord = column / 32;
            int columnBit = 31 - (column % 32);
            if((v[columnWord] & (1 << columnBit)) != 0) {
                for (int rowWord = 0; rowWord < wordsPerColumn; ++rowWord) {
                    ret[rowWord] ^= elements[column * wordsPerColumn + rowWord];
                }
            }
        }
        return ret;
    }

    void setElement(int row, int column) {
        int rowWord = row / 32;
        int rowBit = 31 - (row % 32);
        elements[column * wordsPerColumn + rowWord] |= 1 << rowBit;
    }

    uint32_t getElement(int row, int column) const {
        int rowWord = row / 32;
        int rowBit = 31 - (row % 32);
        return (elements[column * wordsPerColumn + rowWord] & (1 << rowBit)) == 0 ? 0 : 1;
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

int main(int argc, char ** argv) {
    int n = 624;
    int w = 32;
    int r = 31;
    int m = 397;
    uint32_t a = 0x9908b0df;
    int rowOffset = 0;
    int columnOffset = 0;
    int blockStart = w * m;
    SparseMatrix A;
    for(int shiftRow = 0; shiftRow < m * w; ++shiftRow) {
        A.insert(shiftRow, w * (n - m) + shiftRow);
    }
    for(int z = 0; z < n - m; ++z) {
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
         * Handle the upper w-r bits of x_k in multiplication by A
         */
        for(int row = rowOffset; row < rowOffset + w - r; ++row) {
            int col = row;
            A.insert(blockStart + row + 1, col);
        }

        /*
         * Handle the lower r bits of x_(k+1) in multiplication by A
         */
        for(int row = rowOffset + w - r; row < rowOffset + w - 1; ++row) {
            int col = (row + w) % (n*w);
            A.insert(blockStart + row + 1, col);
        }

#if 1
        /*
         * Handle the last column of A
         */
        for(int i = 0; i < w; ++i) {
            int row = rowOffset + i;
            int col = (columnOffset + 2*w - 1) % (n*w);
            // The low order bit of 'a' is in the lower left corner of the matrix, hence a shift of w-i
            if(a & (1 << (w-1-i))) {
                A.insert(blockStart + row, col);
            }
        }

#endif
        /*
         * Handle the sum of x_(k+m)
         */
        for(int i = 0; i < w; ++i) {
            int row = rowOffset + i;
            int col = (columnOffset + m * w + i) % (n*w);
            A.insert(blockStart + row, col);
        }

        rowOffset += w;
        columnOffset += w;
    }

    DenseF2Matrix denseMatrix = A.bake(w);

    MersenneTwister mt19937 = MersenneTwister::createMT19937();
    vector<uint32_t> state = mt19937.getState();

    vector<uint32_t> product = denseMatrix.multiply(state);
    int numIterationsChecked = 1;
    int k = 1000;
    for(int j = 0; j < numIterationsChecked; ++j) {
        for (int i = 0; i < (n - m); ++i) {
            uint32_t r = mt19937.generateUntempered();
            if (product[m + i] != r) {
                cout << "mistake at position " << k << " " << product[m + i] << " " << r << "\n";
            }
            ++k;
        }
        product = denseMatrix.multiply(product);
    }
    cout << "Iteration check done on " << k * (n-m) << " elements" << endl;

    cout << "Num nonzero elements in A " << A.getNumElements() << "\n";

    int numMultiplies = 15;
    timespec t1, t2;
    uint64_t power = 1;
    for(int i = 0; i < numMultiplies; ++i) {
        clock_gettime(CLOCK_REALTIME, &t1);
        denseMatrix = denseMatrix.square();
        power = power + power;
        clock_gettime(CLOCK_REALTIME, &t2);
        double time = (t2.tv_sec * 1E9 + t2.tv_nsec - t1.tv_sec * 1E9 - t1.tv_nsec)/1E9;
        cout << "time for multiply " << i + 1 << ": " << time << " seconds (the matrix is getting denser)\n";
    }
    cout << "power: " << power << "\n";
    cout << "iterate: " << power*(n-m) << "\n";
    return 0;
}