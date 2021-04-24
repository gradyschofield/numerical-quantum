//
// Created by Grady Schofield on 4/24/21.
//

#ifndef ATOM_SPHERICALPOINT_H
#define ATOM_SPHERICALPOINT_H


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

    static double lengthOfSum(SphericalPoint const & p1, SphericalPoint const & p2) {
        double x1 = p1.getR() * cos(p1.getPhi()) * sin(p1.getTheta());
        double y1 = p1.getR() * sin(p1.getPhi()) * sin(p1.getTheta());
        double z1 = p1.getR() * cos(p1.getTheta());

        double x2 = p2.getR() * cos(p2.getPhi()) * sin(p2.getTheta());
        double y2 = p2.getR() * sin(p2.getPhi()) * sin(p2.getTheta());
        double z2 = p2.getR() * cos(p2.getTheta());
        return sqrt( (x1+x2)*(x1+x2) + (y1+y2)*(y1+y2) + (z1+z2)*(z1+z2));
    }
};

#endif //ATOM_SPHERICALPOINT_H
