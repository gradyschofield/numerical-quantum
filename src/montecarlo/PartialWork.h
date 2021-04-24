//
// Created by Grady Schofield on 4/24/21.
//

#ifndef ATOM_PARTIALWORK_H
#define ATOM_PARTIALWORK_H

#include<vector>

using namespace std;

template<typename AccumulatorType>
class alignas(128) PartialWork {
    AccumulatorType inVolume;
    AccumulatorType outVolume;
    AccumulatorType integral;
    double minIntegral = numeric_limits<double>::max();
    double maxIntegral = -numeric_limits<double>::max();
public:
    void addInVolume(double x) {
        inVolume += x;
    }

    void addOutVolume(double x) {
        outVolume += x;
    }

    void addIntegral(double x) {
        integral += x;
        minIntegral = min(minIntegral, x);
        maxIntegral = max(maxIntegral, x);
    }

    double getInVolume() const {
        return inVolume;
    }

    double getOutVolume() const {
        return outVolume;
    }

    double getIntegral() const {
        return integral;
    }

    double getMinIntegral() const {
        return minIntegral;
    }

    double getMaxIntegral() const {
        return maxIntegral;
    }

    static PartialWork accumulate(vector<PartialWork> const & v) {
        PartialWork partialWork;
        for(PartialWork const & p : v) {
            partialWork.addInVolume(p.getInVolume());
            partialWork.addOutVolume(p.getOutVolume());
            partialWork.addIntegral(p.getIntegral());
        }
        return partialWork;
    }
};

#endif //ATOM_PARTIALWORK_H
