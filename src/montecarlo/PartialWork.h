//
// Created by Grady Schofield on 4/24/21.
//

#ifndef ATOM_PARTIALWORK_H
#define ATOM_PARTIALWORK_H

#include<vector>

using namespace std;

class alignas(128) PartialWork {
    double inVolume = 0;
    double outVolume = 0;
    double integral = 0;
public:
    double addInVolume(double x) {
        inVolume += x;
        return inVolume;
    }

    double addOutVolume(double x) {
        outVolume += x;
        return outVolume;
    }

    double addIntegral(double x) {
        integral += x;
        return integral;
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
