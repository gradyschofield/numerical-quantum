//
// Created by grady on 6/30/20.
//

#include<vector>
#include<iostream>
#include<fstream>
#include<complex>
#include<EigenSolver.h>
#include <grid/DiracPhysics.h>

using namespace std;

double const speedOfLight = 274;
double const massOfElectron = 0.5;
double const chargeOfElectron = sqrt(2);

std::vector<double> deriv2{-0.5, 0, 0.5};
std::vector<double> deriv4{1/12.0, -2/3.0, 0, 2/3.0, -1/12.0};
std::vector<double> deriv6{-1/60.0, 3/20.0, -3/4.0, 0, 3/4.0, -3/20.0, 1/60.0};

vector<complex<double>> gamma0(vector<complex<double>> v) {
    return vector<complex<double>>{v[0], v[1], -v[2], -v[3]};
}

vector<complex<double>> gamma1(vector<complex<double>> v) {
    return vector<complex<double>>{v[3], v[2], -v[1], -v[0]};
}

vector<complex<double>> gamma2(vector<complex<double>> v) {
    return vector<complex<double>>{
            complex<double>(v[3].imag(), -v[3].real()),
            complex<double>(-v[2].imag(), v[2].real()),
            complex<double>(-v[1].imag(), v[1].real()),
            complex<double>(v[0].imag(), -v[0].real())
    };
}

vector<complex<double>> gamma3(vector<complex<double>> v) {
    return vector<complex<double>>{v[2], -v[3], -v[0], v[1]};
}

struct Hamiltonian {
    int nx;
    int ny;
    int nz;
    double px, py, pz;
    int ndim;
    double h;
    double hinv;

    Hamiltonian(double boxSize, double h) {
        Hamiltonian::h = h;
        hinv = 1/h;
        nx = ny = nz = boxSize / h;
        px = h/2;
        py = h/2;
        pz = h/2;
        ndim = (1+2*nx) * (1+2*ny) * (1+2*nz);
    }

    int getIndex(int spinorComponent, int i, int j, int k) const {
        int zeroBasedI = i + nx;
        int zeroBasedJ = j + ny;
        int zeroBasedK = k + nz;
        int numYPoints = 1 + 2 * ny;
        int numZPoints = 1 + 2 * nz;
        int offset = zeroBasedI * numYPoints * numZPoints + zeroBasedJ * numZPoints + zeroBasedK;
        return offset + spinorComponent * ndim;
    }

    void getGradients(vector<complex<double>> & gradients, vector<complex<double>> const & v, int i, int j, int k) const {
        int stencilHalfWidth = 3;
        for(int spinorComponent = 0; spinorComponent < 4; ++spinorComponent) {
            complex<double> dx;
            for (int ii = max(i - stencilHalfWidth, -nx); ii <= min(i + stencilHalfWidth, nx); ++ii) {
                int stencilOffset = (ii - i) + stencilHalfWidth;
                int ndx = getIndex(spinorComponent, ii, j, k);
                dx += deriv6[stencilOffset] * v[ndx];
            }
            gradients[spinorComponent*3] = dx * hinv;
            complex<double> dy;
            for (int jj = max(j - stencilHalfWidth, -ny); jj <= min(j + stencilHalfWidth, ny); ++jj) {
                int stencilOffset = (jj - j) + stencilHalfWidth;
                int ndx = getIndex(spinorComponent, i, jj, k);
                dy += deriv6[stencilOffset] * v[ndx];
            }
            gradients[spinorComponent*3 + 1] = dy * hinv;
            complex<double> dz;
            for (int kk = max(k - stencilHalfWidth, -nz); kk <= min(k + stencilHalfWidth, nz); ++kk) {
                int stencilOffset = (kk - k) + stencilHalfWidth;
                int ndx = getIndex(spinorComponent, i, j, kk);
                dz += deriv6[stencilOffset] * v[ndx];
            }
            gradients[spinorComponent*3+2] = dz * hinv;
        }
    }

    void matvec(vector<complex<double>> & ret, vector<complex<double>> const & v) const {
        for(auto & x : ret) {
            x = complex<double>();
        }
        vector<complex<double>> gradients(12);
        int retIdx = 0;
        for (int i = -nx; i <= nx; ++i) {
            for (int j = -ny; j <= ny; ++j) {
                for (int k = -nz; k <= nz; ++k) {
                    getGradients(gradients, v, i, j, k);
                    vector<complex<double>> v1 = gamma0(gamma1({ gradients[0], gradients[3], gradients[6], gradients[9]}));
                    vector<complex<double>> v2 = gamma0(gamma2({ gradients[1], gradients[4], gradients[7], gradients[10]}));
                    vector<complex<double>> v3 = gamma0(gamma3({ gradients[2], gradients[5], gradients[8], gradients[11]}));
                    double x = nx * h;
                    double y = ny * h;
                    double z = nz * h;
                    double r = sqrt((x-px)*(x-px) + (y-py)-(y-py) + (z-pz)*(z-pz));
                    for(int m = 0; m < 4; ++m) {
                        /*
                         *  Kinetic energy
                         */
                        ret[retIdx + ndim*m] -= complex<double>(0, speedOfLight) * (v1[m] + v2[m] + v3[m]);
                        /*
                         *  Rest mass energy
                         */
                        ret[retIdx + ndim*m] += (m < 2 ? 1 : -1) * massOfElectron * speedOfLight * speedOfLight * v[retIdx + ndim*m];
                        /*
                         *  Coulomb potential from nucleus a (px, py, pz)
                         */
                        ret[retIdx + ndim*m] += chargeOfElectron * v[retIdx + ndim*m] / r;
                    }
                    ++retIdx;
                }
            }
        }
    }
};

struct ElectronPositronNorm {
    double electronPart;
    double positronPart;
};

ElectronPositronNorm getElectronPositronParts(vector<complex<double>> const & v) {
    double e = 0;
    double p = 0;
    for(int i = 0; i < v.size()/2; ++i) {
        e += norm(v[i]);
    }
    for(int i = v.size()/2; i < v.size(); ++i) {
        p += norm(v[i]);
    }
    return ElectronPositronNorm{e, p};
}

int main(int argc, char ** argv) {
    /*
    ChebyshevFilter chebyshevFilter = ChebyshevFilter::triangle(200, -1, 1, -0.1, 0.1);
    chebyshevFilter.plot("/tmp/cheb");
    return 0;
     */
    EigenSolver<grid::DiracPhysics> eigenSolver;
    double h = 0.5;
    double boxHalfWidth = 8;
    Hamiltonian hamiltonian(boxHalfWidth, h);
    cout << "Hamiltonian dimension: " << hamiltonian.ndim << "\n";
    auto matvec = [&hamiltonian](vector<complex<double>> & ret, vector<complex<double>> const & v) {
        hamiltonian.matvec(ret, v);
    };

    pair<double, double> bounds = eigenSolver.spectrumBounds(matvec, hamiltonian.ndim);
    cout << "Spectrum bounds: " << bounds.first << " " << bounds.second << "\n";

    vector<vector<complex<double>>> subspace = eigenSolver.generateStepFilteredSubspace(matvec, 250,
                                                                                        5,
                                                                                        hamiltonian.ndim,
                                                                                        bounds,
                                                                                        make_pair(-2000, 2000), 80);
    /*
    for(int i = 0; i < subspace.size(); ++i) {
        for(int j = i; j < subspace.size(); ++j) {
            cout << "dot subspace: " << i << " " << j << " " << dot(subspace[i], subspace[j]) << "\n";
        }
    }
     */
    cout << "Starting Rayleigh-Ritz" << endl;
    EigenSolver<grid::DiracPhysics>::EigenSystem eigenSystem = eigenSolver.rayleighRitz(matvec, subspace);
    for(int i = 0; i < eigenSystem.eigenvalues.size(); ++i) {
        ElectronPositronNorm n = getElectronPositronParts(eigenSystem.eigenvectors[i]);
        cout << eigenSystem.eigenvalues[i] << " " << n.electronPart << " " << n.positronPart << "\n";
    }
    return 0;
}