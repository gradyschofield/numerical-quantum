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

class BigSum {
    vector<double> positiveAccumulators;
    vector<double> negativeAccumulators;
    vector<double> orderLimits;
    int orders;

public:

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
        double lowestLimit = pow(10, round(log10(minBin)));
        for(int i = 0; i < orders; ++i) {
            orderLimits[i] = lowestLimit * pow(10,i);
        }
    }

    void add(double x) {
        int bin = min(orders-1, (int)max(0.0, round(log10(fabs(x)))));
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
};

#endif //ATOM_BIGSUM_H
