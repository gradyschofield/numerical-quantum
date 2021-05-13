//
// Created by Grady Schofield on 5/10/21.
//

#include<iostream>
#include<fstream>
#include<sstream>
#include<unordered_set>
#include<vector>
#include<thread>
#include<tuple>

#include<MersenneTwister.h>
#include<DenseF2Matrix.h>
#include<SparseMatrix.h>

using namespace std;

/*
 * This program implements the fast jump ahead idea in
 * http://www.iro.umontreal.ca/~lecuyer/myftp/papers/jumpf2.pdf .
 * It also implements the slower straightforward way to jump ahead as a reference
 * to compare the fast code for correctness.
 */

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
    int numIterationsChecked = 10;
    int k = 0;
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

    vector<uint32_t> poly = denseMatrix.characteristicPolynomial();

    /*
    int numMultiplies = 19940;
    timespec t1, t2;
    uint64_t power = 1;
    cout << "The matrix gets denser as more multiplications are done, hence the slowdown\n";
    DenseF2Matrix matrixAccum = denseMatrix;
    ofstream traceOfs("traces");
    uint32_t trace = denseMatrix.trace();
    if(trace) {
        traceOfs << "1 " << denseMatrix.trace() << "\n";
    }
    for(int i = 0; i < numMultiplies; ++i) {
        clock_gettime(CLOCK_REALTIME, &t1);
        matrixAccum = denseMatrix.multiply(matrixAccum);
        clock_gettime(CLOCK_REALTIME, &t2);
        double time = (t2.tv_sec * 1E9 + t2.tv_nsec - t1.tv_sec * 1E9 - t1.tv_nsec)/1E9;
        cout << "time for multiply " << i + 1 << ": " << time << " seconds.\n";
        trace = matrixAccum.trace();
        if(trace) {
            traceOfs << i + 2 << " " << matrixAccum.trace() << "\n";
            traceOfs.flush();
        }

    }
    cout << "power: " << power << "\n";
    cout << "iterate: " << power*(n-m) << "\n";
     */
    return 0;
}