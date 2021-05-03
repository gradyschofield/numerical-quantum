//
// Created by Grady Schofield on 4/24/21.
//

#ifndef ATOM_SPHERICALPOINT_H
#define ATOM_SPHERICALPOINT_H

#include<cmath>

#include<EuclideanPoint.h>
#include<MathUtil.h>

class SphericalPoint {
    double r, theta, phi;
public:
    SphericalPoint(double r, double theta, double phi) : r(r), theta(theta), phi(phi){}

    double getR() const {
        return r;
    }

    double getTheta() const {
        return theta;
    }

    double getPhi() const {
        return phi;
    }

    double volumeElementAt() const {
        return getR() * getR() * sin(getTheta());
    }

    EuclideanPoint toEuclideanPoint() const {
        double x = getR() * cos(getPhi()) * sin(getTheta());
        double y = getR() * sin(getPhi()) * sin(getTheta());
        double z = getR() * cos(getTheta());
        return EuclideanPoint(x, y, z);
    }

    double dot(SphericalPoint const & p) const {
        return toEuclideanPoint().dot(p.toEuclideanPoint());
    }

    SphericalPoint operator+(SphericalPoint const & p) const {
        return fromEuclidean(toEuclideanPoint() + p.toEuclideanPoint());
    }

    static SphericalPoint fromEuclidean(EuclideanPoint const & p) {
        double r = sqrt(square(p.getX()) + square(p.getY()) + square(p.getZ()));
        if(r == 0) {
            return SphericalPoint(0, 0, 0);
        }
        double theta = acos(p.getZ() / r);
        if(theta == 0 || theta == M_PI) {
            return SphericalPoint(r, theta, 0);
        }
        double phi = atan2(p.getY(), p.getX());
        return SphericalPoint(r, theta, phi);
    }

    static double lengthOfSum(SphericalPoint const & p1, SphericalPoint const & p2) {
        return (p1.toEuclideanPoint() + p2.toEuclideanPoint()).norm();
    }

    static double lengthOfSum(vector<SphericalPoint> const & pointList) {
        EuclideanPoint e;
        for(auto & p : pointList) {
            e += p.toEuclideanPoint();
        }
        return e.norm();
    }
};

ostream & operator<<(ostream & ofs, SphericalPoint const & p) {
    ofs << "(" << p.getR() << ", " << p.getTheta() << ", " << p.getPhi() << ")";
    return ofs;
}

#endif //ATOM_SPHERICALPOINT_H
