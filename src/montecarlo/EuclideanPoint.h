//
// Created by Grady Schofield on 4/25/21.
//

#ifndef ATOM_EUCLIDEANPOINT_H
#define ATOM_EUCLIDEANPOINT_H

class EuclideanPoint {
    double x, y, z;
public:

    EuclideanPoint()
        : x(0), y(0), z(0)
    {
    }

    EuclideanPoint(double x, double y, double z) : x(x), y(y), z(z) {}

    double getX() const {
        return x;
    }

    double getY() const {
        return y;
    }

    double getZ() const {
        return z;
    }

    double dot(EuclideanPoint const & e) const {
        return x*e.x + y*e.y + z*e.z;
    }

    double norm() const {
        return sqrt(x*x + y*y + z*z);
    }

    EuclideanPoint operator+(EuclideanPoint p) const {
        return EuclideanPoint(x+p.x, y+p.y, z+p.z);
    }

    void operator+=(EuclideanPoint const & p) {
        x += p.x;
        y += p.y;
        z += p.z;
    }
};

#endif //ATOM_EUCLIDEANPOINT_H
