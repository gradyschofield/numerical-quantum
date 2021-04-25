//
// Created by grady on 4/24/21.
//

#ifndef ATOM_BIGSUM_H
#define ATOM_BIGSUM_H

#include<vector>
#include<cmath>
#include<algorithm>
#include<iomanip>

using namespace std;

/*
 * The point of this class is to avoid cancellation and truncation errors in sums that may be
 * composed of trillions of terms and where the relative sizes of the terms is not known.
 * It's meant to be used in Monte Carlo integration.
 *
 * This class represents a sum by keeping partial results in several bins, each of which stores a
 * number that doesn't exceed a certain order of magnitude. When a number is added to the sum,
 * it is accumulated into the appropriate bin based on its size.  When the value in a bin goes
 * over the order of magnitude limit for that bin, the bin value is added to the bin for next largest
 * order of magnitude and the original bin is set to zero.  This "carry" operation can cause a cascade
 * of carries all the  way up to the highest bin.
 *
 * The point of this is to avoid adding small numbers to large ones until the small ones have accumulated
 * into a more substantial quantity.  Negative and positive terms are accumulated separately to
 * avoid cancellation.
 */

class BigSum {
    vector<double> positiveAccumulators;
    vector<double> negativeAccumulators;
    vector<double> orderLimits;
    int orders;
    int binOffset;

    void carryAccumulators(int startBin, bool positive) {
        if(positive) {
            for (int bin = startBin; bin < orders - 1; ++bin) {
                if (positiveAccumulators[bin] > orderLimits[bin]) {
                    positiveAccumulators[bin + 1] += positiveAccumulators[bin];
                    positiveAccumulators[bin] = 0;
                } else {
                    break;
                }
            }
        } else {
            for (int bin = startBin; bin < orders - 1; ++bin) {
                if (negativeAccumulators[bin] < -orderLimits[bin]) {
                    negativeAccumulators[bin + 1] += negativeAccumulators[bin];
                    negativeAccumulators[bin] = 0;
                } else {
                    break;
                }
            }
        }
    }

public:

    BigSum()
            : BigSum(25, 1E-5)
    {
    }

    BigSum(int orders, double minBin)
            : positiveAccumulators(orders),
              negativeAccumulators(orders),
              orderLimits(orders),
              orders(orders)
    {
        binOffset = round(log10(minBin));
        double lowestLimit = pow(10, binOffset);
        for(int i = 0; i < orders; ++i) {
            orderLimits[i] = lowestLimit * pow(10,i);
        }
    }

    operator double() const {
        double positiveSum = 0;
        double negativeSum = 0;
        for(int i = 0; i < orders; ++i) {
            positiveSum += positiveAccumulators[i];
        }
        for(int i = 0; i < orders; ++i) {
            negativeSum += negativeAccumulators[i];
        }
        return positiveSum + negativeSum;
    }

    void add(double x) {
        int bin = min(orders-1, (int)max(0.0, round(log10(fabs(x))) - binOffset));
        if(x > 0) {
            positiveAccumulators[bin] += x;
            if (positiveAccumulators[bin] > orderLimits[bin]) {
                carryAccumulators(bin, true);
            }
        } else {
            negativeAccumulators[bin] += x;
            if (negativeAccumulators[bin] < -orderLimits[bin]) {
                carryAccumulators(bin, false);
            }
        }
    }

    void operator+=(double x) {
        add(x);
    }

    void operator-=(double x) {
        add(-x);
    }

    void print() {
        cout << "orders: " << orders << "\n";
        double sum = 0;
        for(int i = 0; i < orders; ++i) {
            cout << positiveAccumulators[i] << " ";
            sum += positiveAccumulators[i];
        }
        cout << "\n";
        for(int i = 0; i < orders; ++i) {
            cout << negativeAccumulators[i] << " ";
            sum += negativeAccumulators[i];
        }
        cout << "\n";
        for(int i = 0; i < orders; ++i) {
            cout << orderLimits[i] << " ";
        }
        cout << "\n";
        cout << "sum: " << setprecision(16) << sum << setprecision(5) << "\n";
    }

};

#endif //ATOM_BIGSUM_H
