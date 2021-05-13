//
// Created by Grady Schofield on 5/12/21.
//

#ifndef ATOM_DENSEF2MATRIX_H
#define ATOM_DENSEF2MATRIX_H

#include<vector>
#include<thread>
#include<unordered_set>

using namespace std;

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

    DenseF2Matrix multiply(DenseF2Matrix const & m) const {
        DenseF2Matrix ret(dim);
        vector<thread> threads;
        int numThreads = thread::hardware_concurrency();
        int stop = 0;
        for(int i = 0; i < numThreads; ++i) {
            int start = stop;
            stop = start + (dim / numThreads) + (i < dim % numThreads ? 1 : 0);
            threads.emplace_back([start, stop, &m, &ret, this]() {
                for (int column = start; column < stop; ++column) {
                    ret.setColumn(column, multiply(m.getColumn(column)));
                }
            });
        }
        for(thread & t : threads) {
            t.join();
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

    uint32_t trace() const {
        uint32_t accum = 0;
        for(int column = 0; column < dim; ++column) {
            int row = column;
            int rowWord = row / 32;
            int rowBit = 31 - row % 32;
            accum ^= (elements[column * wordsPerColumn + rowWord] & (1 << rowBit)) == 0 ? 0 : 1;
        }
        return accum;
    }

    static vector<uint32_t> polyProduct(bool plusOne, vector<uint32_t> polynomial) {
        uint32_t carry = polynomial[0] & 0x80000000;
        polynomial[0] = polynomial[0] << 1 | (plusOne ? polynomial[0] : 0);
        for(int i = 1; i < polynomial.size(); ++i) {
            uint32_t nextCarry = polynomial[i] & 0x80000000;
            polynomial[i] = polynomial[i] << 1 | (plusOne ? polynomial[i] : 0);
            if(carry) {
                polynomial[i] |= 0x1;
            }
            carry = nextCarry;
        }
        return polynomial;
    }

    static void polySum(vector<uint32_t> & totalPoly, vector<uint32_t> const & polynomial) {
        for(int i = 0; i < polynomial.size(); ++i) {
            totalPoly[i] ^= polynomial[i];
        }
    }

    vector<uint32_t> characteristicPolynomialRecurse(uint32_t & sums, int row, unordered_set<int> & skipRows, vector<uint32_t> prefactor) const {
        static uint32_t lastSum = 0;
        vector<uint32_t> totalPoly(wordsPerColumn);
        if(sums > 0 && sums != lastSum) {
            cout << "Sums: " << sums << "\n";
            lastSum = sums;
        }
        for(int column = 0; column < dim; ++column) {
            if(!skipRows.contains(column) &&
               ( getElement(row, column) == 1 || row == column)) {
                skipRows.insert(column);
                if(row == column) {
                    bool plusOne = 1 == getElement(row, column);
                    polySum(totalPoly, characteristicPolynomialRecurse(sums, row + 1, skipRows, polyProduct(plusOne, prefactor)));
                    sums += 1;
                } else {
                    polySum(totalPoly, characteristicPolynomialRecurse(sums, row + 1, skipRows, prefactor));
                    sums += 1;
                }
                skipRows.erase(column);
            }
        }
        return totalPoly;
    }

    vector<uint32_t> characteristicPolynomial() const {
        vector<uint32_t> prefactor(wordsPerColumn);
        prefactor[0] = 1;
        unordered_set<int> skipRows;
        uint32_t sums = 0;
        return characteristicPolynomialRecurse(sums, 0, skipRows, prefactor);
    }
};

#endif //ATOM_DENSEF2MATRIX_H
