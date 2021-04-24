//
// Created by Grady Schofield on 4/24/21.
//

#ifndef ATOM_GENERATOR_H
#define ATOM_GENERATOR_H

#include<random>

#include<Bounds.h>
#include<SphericalPoint.h>

using namespace std;

class Generator {
    mt19937_64 generator;
    uniform_real_distribution<> unifx;
    uniform_real_distribution<> unify;
    uniform_real_distribution<> unifz;
public:
    Generator(Bounds const & bounds) {
        random_device r;
        seed_seq seed{r(), r(), r(), r(), r(), r(), r(), r()};
        generator = mt19937_64(seed);
        unifx = uniform_real_distribution<>(bounds.getMinX(), bounds.getMaxX());
        unify = uniform_real_distribution<>(bounds.getMinY(), bounds.getMaxY());
        unifz = uniform_real_distribution<>(bounds.getMinZ(), bounds.getMaxZ());
    }

    SphericalPoint generate() {
        return SphericalPoint(unifx(generator), unify(generator), unifz(generator));
    }
};

#endif //ATOM_GENERATOR_H
