//
// Created by grady on 7/19/20.
//

#ifndef ATOM_DIRACPHYSICS_H
#define ATOM_DIRACPHYSICS_H

#include<vector>
#include<cstdlib>
#include<cmath>
#include<complex>

using namespace std;

namespace grid {
    class DiracPhysics {
    public:
        typedef complex<double> scalarType;

        static vector<complex<double>> randomVector(int ndim) {
            vector<complex<double>> ret(4 * ndim);
            double rinv = 1.0/RAND_MAX;
            double n = 0;
            for(int i = 0; i < ndim; ++i) {
                vector<complex<double>> v(4);
                for(int j = 0; j < 4; ++j) {
                    v[j].real(rand()*rinv);
                    v[j].imag(rand()*rinv);
                    ret[i*4 + j] = v[j];
                }
                n += std::norm(v[0]) + std::norm(v[1]) + std::norm(v[2]) + std::norm(v[3]);
            }
            complex<double> ninv = 1/sqrt(n);
            for(complex<double> & x : ret) {
                x *= ninv;
            }
            return ret;
        }

        static double norm(vector<complex<double>> const & v) {
            double n = 0;
            for(int i = 0; i < v.size(); i += 4) {
                vector<complex<double>> t(4);
                for(int j = 0; j < 4; ++j) {
                    t[j] = v[i + j];
                }
                n += std::norm(t[0]) + std::norm(t[1]) + std::norm(t[2]) + std::norm(t[3]);
            }
            return sqrt(n);
        }

        static complex<double> dot(vector<complex<double>> const & v1, vector<complex<double>> const & v2) {
            complex<double> d;
            for(int i = 0; i < v2.size(); ++i) {
                d += conj(v1[i]) * v2[i];
            }
            return d;
        }

        static void saxpy(vector<complex<double>> & x, complex<double> a, vector<complex<double>> const & y) {
            for(int i = 0; i < x.size(); ++i) {
                x[i] += a * y[i];
            }
        }

        static double normalize(vector<complex<double>> & v) {
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

#endif //ATOM_DIRACPHYSICS_H
