//
// Created by grady on 7/21/20.
//

#ifndef ATOM_REALSCALARANTISYMMETRICPHYSICS_H
#define ATOM_REALSCALARANTISYMMETRICPHYSICS_H

#include<grid/RealScalarPhysics.h>

namespace grid {
    class RealScalarAntisymmetricPhysics : public RealScalarPhysics {
    public:
        static double & elem(vector<double> & v, int i, int j, int ndim) {
            return v[i + j * ndim];
        }

        static double elem(vector<double> const & v, int i, int j, int ndim) {
            return v[i + j * ndim];
        }

        static double readElem(vector<double> const & v, int i, int j, int ndim) {
            if(i == 0) {
                return 0;
            } else if(i < j) {
                return v[i + j * ndim];
            } else {
                return -v[j + i * ndim];
            }
        }

        static void antiSymmetrize(vector<double> & v, int ndim) {
            return;
            for(int i = 0; i < ndim; ++i) {
                for(int j = i; j < ndim; ++j) {
                    if(i == j) {
                        elem(v, i, j, ndim) = 0;
                    } else {
                        elem(v, j, i, ndim) = -elem(v, i, j, ndim);
                    }
                }
            }
        }

        static vector<double> randomVector(int ndim) {
            vector<double> ret(ndim*ndim);
            double rinv = 1.0/RAND_MAX;
            double n = 0;
            for(int i = 0; i < ndim; ++i) {
                for(int j = i; j < ndim; ++j) {
                    if(i == j) {
                        elem(ret, i, j, ndim) = 0;
                    } else {
                        double t = rand() * rinv;
                        elem(ret, i, j, ndim) = t;
                        elem(ret, j, i, ndim) = -t;
                        n += 2 * t * t;
                    }
                }
            }
            double ninv = 1/sqrt(n);
            for(double & x : ret) {
                x *= ninv;
            }
            return ret;
        }
    };
}

#endif //ATOM_REALSCALARANTISYMMETRICPHYSICS_H
