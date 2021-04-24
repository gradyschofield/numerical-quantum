//
// Created by Grady Schofield on 4/24/21.
//

#ifndef ATOM_BOUNDS_H
#define ATOM_BOUNDS_H

#include<limits>
#include<cmath>

using namespace std;

class alignas(128) Bounds {
    double minx = numeric_limits<double>::max();
    double maxx = -numeric_limits<double>::max();
    double miny = numeric_limits<double>::max();
    double maxy = -numeric_limits<double>::max();
    double minz = numeric_limits<double>::max();
    double maxz = -numeric_limits<double>::max();

public:
    Bounds(){
    }

    Bounds(double minx, double maxx,
           double miny, double maxy,
           double minz, double maxz)
            : minx(minx), maxx(maxx), miny(miny), maxy(maxy), minz(minz), maxz(maxz)
    {
    }

    double getMinX() const {
        return minx;
    }

    double getMaxX() const {
        return maxx;
    }

    double getMinY() const {
        return miny;
    }

    double getMaxY() const {
        return maxy;
    }

    double getMinZ() const {
        return minz;
    }

    double getMaxZ() const {
        return maxz;
    }

    double getVolume() const {
        return (pow(maxx,3) - pow(minx, 3))/3 * (cos(miny) - cos(maxy)) * (maxz - minz);
    }
};

#endif //ATOM_BOUNDS_H
