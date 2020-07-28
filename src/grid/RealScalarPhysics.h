//
// Created by grady on 7/19/20.
//

#ifndef ATOM_REALSCALARPHYSICS_H
#define ATOM_REALSCALARPHYSICS_H

#include<vector>
#include<cstdlib>
#include<cmath>

using namespace std;

namespace grid {
    /*
     *  Must define saxpy, dot, norm, normalize, matvec, randomVector
     */
    class RealScalarPhysics {
    public:
        typedef double scalarType;

        static void antiSymmetrize(vector<double> & v, int ndim) {
        }

        static vector<double> randomVector(int ndim) {
            vector<double> ret(ndim);
            double rinv = 1.0/RAND_MAX;
            double n = 0;
            for(int i = 0; i < ndim; ++i) {
                ret[i] = rand()*rinv;
                n += ret[i] * ret[i];
            }
            double ninv = 1/sqrt(n);
            for(double & x : ret) {
                x *= ninv;
            }
            return ret;
        }

        static double norm(vector<double> const & v) {
            double n = 0;
            for(int i = 0; i < v.size(); ++i) {
                n += v[i] * v[i];
            }
            return sqrt(n);
        }

        static double dot(vector<double> const & v1, vector<double> const & v2) {
            double d = 0;
            for(int i = 0; i < v2.size(); ++i) {
                d += v1[i] * v2[i];
            }
            return d;
        }

        static void saxpy(vector<double> & x, double a, vector<double> const & y) {
            for(int i = 0; i < x.size(); ++i) {
                x[i] += a * y[i];
            }
        }

        static double normalize(vector<double> & v) {
            double n = norm(v);
            if(n == 0) {
                return 0;
            }
            double ninv = 1.0/n;
            for(auto & x : v) {
                x *= ninv;
            }
            return n;
        }
    };
}
#endif //ATOM_REALSCALARPHYSICS_H
