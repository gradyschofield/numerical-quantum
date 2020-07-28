//
// Created by grady on 7/19/20.
//

#ifndef ATOM_EIGENSOLVER_H
#define ATOM_EIGENSOLVER_H

#include<vector>
#include<iostream>
#include<lapacke.h>
#include "ChebyshevFilter.h"

using namespace std;

template<typename Physics>
class EigenSolver{
public:
    void orth(vector<vector<typename Physics::scalarType>> const & subspace, vector<typename Physics::scalarType> &v) const {
        Physics::normalize(v);
        int passes = 0;
        while(true) {
            for (int i = 0; i < subspace.size(); ++i) {
                typename Physics::scalarType d = Physics::dot(subspace[i], v);
                Physics::saxpy(v, -d, subspace[i]);
                Physics::antiSymmetrize(v, sqrt(v.size()));
            }
            ++passes;
            double endNorm = Physics::normalize(v);
            if(endNorm > 1/2 && passes > 1) {
                break;
            }
        }
        /*
        for (int i = 0; i < subspace.size(); ++i) {
            typename Physics::scalarType d = Physics::dot(subspace[i], v);
            cout << "check orth: " << d << "\n";
        }
         */
    }

    template<typename Matvec>
    pair<double, double> spectrumBounds(Matvec && matvec, int ndim, int subspaceSize = 50) const {
        subspaceSize = std::min(ndim, subspaceSize);
        vector<vector<typename Physics::scalarType>> subspace = generateSubspace(matvec, subspaceSize, ndim);
        EigenSystem eigenSystem = rayleighRitz(matvec, subspace);
        return make_pair(eigenSystem.eigenvalues.front(), eigenSystem.eigenvalues.back());
    }

    template<typename Matvec>
    vector<vector<typename Physics::scalarType>> generateSubspace(Matvec && matvec, int subspaceSize, int ndim) const {
        vector<vector<typename Physics::scalarType>> subspace;
        vector<typename Physics::scalarType> v = Physics::randomVector(ndim);
        for(int i = 0; i < subspaceSize; ++i) {
            vector<typename Physics::scalarType> ret(v.size());
            matvec(ret, v);
            orth(subspace, ret);
            v = ret;
            subspace.emplace_back(move(ret));
        }
        return subspace;
    }

    template<typename Matvec>
    void scaledMatvec(Matvec && matvec, double shift, double scale, vector<typename Physics::scalarType> & ret, vector<typename Physics::scalarType> const & v) const {
        matvec(ret, v);
        for(int i = 0; i < ret.size(); ++i) {
            ret[i] = scale * (ret[i] - shift*v[i]);
        }
        Physics::antiSymmetrize(ret, sqrt(ret.size()));
    }

    template<typename Matvec>
    void applyFilter(Matvec && matvec, ChebyshevFilter const & chebyshevFilter, pair<double, double> bounds, vector<typename Physics::scalarType> & ret, vector<typename Physics::scalarType> const & v) const {
        double shift = (bounds.first + bounds.second) / 2;
        double scale = 2.0 / (1.01*(bounds.second - bounds.first)); // expand the size of the interval by 1% for safety
        vector<typename Physics::scalarType> iterate1(v.size());
        vector<typename Physics::scalarType> iterate2(v.size());
        vector<typename Physics::scalarType> iterate3(v.size());
        vector<typename Physics::scalarType> t(v.size());
        double factor = chebyshevFilter.coef[0] / sqrt(M_PI);
        for(int i = 0; i < ret.size(); ++i) {
            iterate1[i] = v[i];
            ret[i] = factor * iterate1[i];
        }
        scaledMatvec(matvec, shift, scale, iterate2, v);
        factor = chebyshevFilter.coef[1] / sqrt(M_PI_2);
        for(int i = 0; i < ret.size(); ++i) {
            ret[i] += factor * iterate2[i];
        }
        for(int j = 2; j < chebyshevFilter.coef.size(); ++j) {
            scaledMatvec(matvec, shift, scale, t, iterate2);
            factor = chebyshevFilter.coef[j] / sqrt(M_PI_2);
            for(int i = 0; i < ret.size(); ++i) {
                iterate3[i] = 2.0 * t[i] - iterate1[i];
                ret[i] += iterate3[i] * factor;
            }
            swap(iterate1, iterate2);
            swap(iterate2, iterate3);
        }
        Physics::antiSymmetrize(ret, sqrt(ret.size()));
    }

    template<typename Matvec>
    vector<vector<typename Physics::scalarType>> generateTriangleFilteredSubspace(Matvec && matvec, int subspaceSize, int blockSize,
                                                                                  int ndim,
                                                                                  pair<double, double> const & spectrumBounds,
                                                                                  pair<double, double> stepBounds, int filterOrder) const {
        ChebyshevFilter chebyshevFilter = ChebyshevFilter::triangle(filterOrder, spectrumBounds.first, spectrumBounds.second,
                                                                    stepBounds.first, stepBounds.second);
        chebyshevFilter.plot("/tmp/cheb");
        return generateFilteredSubspace(matvec, subspaceSize, blockSize, ndim, chebyshevFilter);
    }

    template<typename Matvec>
    vector<vector<typename Physics::scalarType>> generateStepFilteredSubspace(Matvec && matvec, int subspaceSize, int blockSize,
                                                                              int ndim,
                                                                              pair<double, double> const & spectrumBounds,
                                                                              pair<double, double> stepBounds, int filterOrder) const {
        ChebyshevFilter chebyshevFilter = ChebyshevFilter::step(filterOrder, spectrumBounds.first, spectrumBounds.second,
                                                                    stepBounds.first, stepBounds.second);
        chebyshevFilter.plot("/tmp/cheb");
        return generateFilteredSubspace(matvec, subspaceSize, blockSize, ndim, chebyshevFilter);
    }

    template<typename Matvec>
    vector<vector<typename Physics::scalarType>> generateFilteredSubspace(Matvec && matvec, int subspaceSize, int blockSize, int ndim, ChebyshevFilter const & chebyshevFilter) const {
        vector<vector<typename Physics::scalarType>> subspace;
        for(int i = 0; i < blockSize; ++i) {
            vector<typename Physics::scalarType> t = Physics::randomVector(ndim);
            vector<typename Physics::scalarType> ret(t.size());
            applyFilter(matvec, chebyshevFilter, chebyshevFilter.spectrumBounds, ret, t);
            orth(subspace, ret);
            subspace.push_back(move(ret));
        }
        for(int i = blockSize; i < subspaceSize; i += blockSize) {
            int startIndex = i - blockSize;
            for(int j = startIndex; j < startIndex + blockSize; ++j) {
                vector<typename Physics::scalarType> v(subspace[j]);
                vector<typename Physics::scalarType> ret(v.size());
                applyFilter(matvec, chebyshevFilter, chebyshevFilter.spectrumBounds, ret, v);
                //cout << "filtered iterate " << i + 1 << "\n";
                orth(subspace, ret);
                subspace.emplace_back(move(ret));
            }
        }
        return subspace;
    }

    struct EigenSystem {
        vector<double> eigenvalues;
        vector<vector<typename Physics::scalarType>> eigenvectors;
    };

    vector<double> hermitianDenseEigenSystem(vector<typename Physics::scalarType> & v) const {
        int dim = (int)sqrt((double)v.size());
        char jobz = 'V';
        char uplo = 'L';
        vector<double> eigenvalues(dim);
        if(is_same<typename Physics::scalarType, complex<double>>::value) {
            LAPACKE_zheev(LAPACK_COL_MAJOR, jobz, uplo, dim, (double _Complex *) v.data(), dim, eigenvalues.data());
        } else if(is_same<typename Physics::scalarType, double>::value) {
            LAPACKE_dsyev(LAPACK_COL_MAJOR, jobz, uplo, dim, (double*)v.data(), dim, eigenvalues.data());
        } else {
            throw runtime_error("Only setup to deal with double and complex<double> in hermitianDenseEigenSystem");
        }
        return eigenvalues;
    }

    template<typename Matvec>
    EigenSystem rayleighRitz(Matvec && matvec, vector<vector<typename Physics::scalarType>> const & subspace) const {
        if(subspace.empty()) {
            throw std::runtime_error("Cannot do Rayleigh-Ritz on empty subspace");
        }
        int ndim = subspace.front().size();
        vector<typename Physics::scalarType> dotMatrix(subspace.size() * subspace.size());
        for(int i = 0; i < subspace.size(); ++i) {
            vector<typename Physics::scalarType> const & v = subspace[i];
            vector<typename Physics::scalarType> av(v.size());
            matvec(av, v);
            for(int j = i; j < subspace.size(); ++j) {
                dotMatrix[j + i * subspace.size()] = Physics::dot(subspace[j], av);
            }
        }
        //cout << "Starting dense solver" << endl;
        EigenSystem ret;
        ret.eigenvalues = hermitianDenseEigenSystem(dotMatrix);

        for(int i = 0; i < subspace.size(); ++i) {
            vector<typename Physics::scalarType> eigv(ndim);
            for(int j = 0; j < subspace.size(); ++j) {
                typename Physics::scalarType c = dotMatrix[j + i * subspace.size()];
                vector<typename Physics::scalarType> const & v = subspace[j];
                for(int k = 0; k < ndim; ++k) {
                    eigv[k] += c * v[k];
                }
            }
            //cout << "eigenvector norm: " << (i+1) << " "  << Physics::norm(eigv) << "\n";
            ret.eigenvectors.emplace_back(move(eigv));
        }
        return ret;
    }

};
#endif //ATOM_EIGENSOLVER_H
