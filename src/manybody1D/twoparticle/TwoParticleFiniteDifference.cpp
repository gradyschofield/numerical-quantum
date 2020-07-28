//
// Created by grady on 7/21/20.
//

#include<grid/RealScalarPhysics.h>
#include<EigenSolver.h>
#include<Constants.h>
#include <grid/RealScalarAntisymmetricPhysics.h>
#include <algorithm>
#include <Energy.h>

vector<double> deriv2_2{1, -2, 1};
vector<double> deriv2_4{-1.0/12, 4.0/3, -5.0/2, 4.0/3, -1.0/12};

double ion1R = -0.4;
double ion2R = -ion1R;

class Hamiltonian {
public:
    int numGridPointsPerDimension;
    int totalGridPoints;
    double h;
    vector<double> deriv2coef;
    int deriv2Length;
    int deriv2BackOffset;
    vector<double> gridPoints;

    Hamiltonian(double gridSpacing, double bounds) {
        h = gridSpacing;
        int ndim = 2 * bounds / gridSpacing;
        for(int i = ndim/2-1; i >= 0; --i) {
            gridPoints.push_back(-h/2 - i * h);
        }
        for(int i = 0; i < ndim/2; ++i) {
            gridPoints.push_back(h/2 + i * h);
        }
        numGridPointsPerDimension = gridPoints.size();
        totalGridPoints = numGridPointsPerDimension * numGridPointsPerDimension;
        deriv2coef = deriv2_2;
        for(auto & coef : deriv2coef) {
            coef *= -1/(2* constants::massOfElectron * h * h);
        }
        for(auto & coef : deriv2coef) {
            cout << "coef: " << coef << "\n";
        }
        deriv2Length = deriv2coef.size();
        deriv2BackOffset = (deriv2Length-1)/2;
        ofstream ofs("/tmp/pot");
        for(double x : gridPoints) {
            ofs << x << " " << potential(x, 1, 0) << "\n";
        }
    }

    double potential(double x, double N, double R) {
        double t = fabs(x - R);
        if(t < 0.5) {
            return constants::chargeOfElectron*N*(4*t*t - 3);
        } else {
            return -constants::chargeOfElectron*N/t;
        }
    }

    /*
     *  This is the force for like charges.  Use -force when integrating election-ion contributions.
     */
    double force(double x, double N, double R) {
        double t = fabs(x - R);
        double sign = (x-R) < 0 ? -1 : 1;
        if(t < 0.5) {
            return -sign*constants::chargeOfElectron*N*(8*t);
        } else {
            return -sign*constants::chargeOfElectron*N/(t*t);
        }
    }

    double ionForce(double x, double N, double R) {
        double t = fabs(x - R);
        double sign = (x-R) < 0 ? -1 : 1;
        return -sign*constants::chargeOfElectron*N/(t*t);
    }

    void matvecKinetic(vector<double> & ret, vector<double> const & v) {
        for(int i = 0; i < numGridPointsPerDimension; ++i) {
            for(int j = 0; j < numGridPointsPerDimension; ++j) {
                for (int k = std::max(0, i - deriv2BackOffset); k <= std::min(numGridPointsPerDimension - 1, i + deriv2BackOffset); ++k) {
                    grid::RealScalarAntisymmetricPhysics::elem(ret, i, j, numGridPointsPerDimension) +=
                            deriv2coef[k - i + deriv2BackOffset] * grid::RealScalarAntisymmetricPhysics::elem(v, k, j, numGridPointsPerDimension);
                }
                for (int k = std::max(0, j - deriv2BackOffset); k <= std::min(numGridPointsPerDimension - 1, j + deriv2BackOffset); ++k) {
                    grid::RealScalarAntisymmetricPhysics::elem(ret, i, j, numGridPointsPerDimension) +=
                            deriv2coef[k - j + deriv2BackOffset] * grid::RealScalarAntisymmetricPhysics::elem(v, i, k, numGridPointsPerDimension);
                }
            }
        }
    }

    /*
     * Ion-Electron interaction
     */
    void matvecIonElectron(vector<double> & ret, vector<double> const & v) {
        for(int i = 0; i < numGridPointsPerDimension; ++i) {
            for (int j = 0; j < numGridPointsPerDimension; ++j) {
                grid::RealScalarAntisymmetricPhysics::elem(ret, i, j, numGridPointsPerDimension) +=
                        grid::RealScalarAntisymmetricPhysics::elem(v, i, j, numGridPointsPerDimension) *
                        (potential(gridPoints[i], 1, ion2R) + potential(gridPoints[i], 1, ion1R) +
                         potential(gridPoints[j], 1, ion2R) + potential(gridPoints[j], 1, ion1R));
            }
        }
    }

    /*
     * Electron-electron interaction
     */
    void matvecElectronElectron(vector<double> & ret, vector<double> const & v) {
        for(int i = 0; i < numGridPointsPerDimension; ++i) {
            for (int j = 0; j < numGridPointsPerDimension; ++j) {
                grid::RealScalarAntisymmetricPhysics::elem(ret, i, j, numGridPointsPerDimension) +=
                        grid::RealScalarAntisymmetricPhysics::elem(v, i, j, numGridPointsPerDimension) *
                        (-potential(gridPoints[i], 1, gridPoints[j]));
            }
        }
    }

    void matvec(vector<double> & ret, vector<double> const & v) {
        for(auto & x : ret) {
            x = 0;
        }
        /*
         * Kinetic energy
         */
        matvecKinetic(ret, v);

        /*
         * Potential Energy
         */
        matvecIonElectron(ret, v);
        matvecElectronElectron(ret, v);
    }

    template<typename Matvec>
    Energy getEnergies(Matvec && matvec, EigenSolver<grid::RealScalarAntisymmetricPhysics>::EigenSystem const & eigenSystem) {
        vector<double> ret(totalGridPoints);
        matvecKinetic(ret, eigenSystem.eigenvectors[0]);
        double kineticEnergy = grid::RealScalarAntisymmetricPhysics::dot(eigenSystem.eigenvectors[0], ret);

        for(auto & x : ret) x = 0;
        matvecIonElectron(ret, eigenSystem.eigenvectors[0]);
        double ionElectronEnergy = grid::RealScalarAntisymmetricPhysics::dot(eigenSystem.eigenvectors[0], ret);

        for(auto & x : ret) x = 0;
        matvecElectronElectron(ret, eigenSystem.eigenvectors[0]);
        double electronElectronEnergy = grid::RealScalarAntisymmetricPhysics::dot(eigenSystem.eigenvectors[0], ret);

        double ionIonEnergy = -potential(ion1R, 1, ion2R);
        return Energy(kineticEnergy, ionElectronEnergy, electronElectronEnergy, 0, ionIonEnergy);
    }
};

double symmetricMeasure(vector<double> const & v, int numGridPointsPerDimension) {
    double measure = 0;
    for(int i = 0; i < numGridPointsPerDimension; ++i) {
        for(int j = i+1; j < numGridPointsPerDimension; ++j) {
            double t1 = grid::RealScalarAntisymmetricPhysics::elem(v, i, j, numGridPointsPerDimension);
            double t2 = grid::RealScalarAntisymmetricPhysics::elem(v, j, i, numGridPointsPerDimension);
            measure += pow(t1-t2, 2);
        }
    }
    return sqrt(measure);
}

double antiSymmetricMeasure(vector<double> const & v, int numGridPointsPerDimension) {
    double measure = 0;
    for(int i = 0; i < numGridPointsPerDimension; ++i) {
        for(int j = i; j < numGridPointsPerDimension; ++j) {
            double t1 = grid::RealScalarAntisymmetricPhysics::elem(v, i, j, numGridPointsPerDimension);
            if(i == j) {
                measure += pow(t1, 2);
            } else {
                double t2 = -grid::RealScalarAntisymmetricPhysics::elem(v, j, i, numGridPointsPerDimension);
                measure += pow(t1 - t2, 2);
            }
        }
    }
    return sqrt(measure);
}

int main(int argc, char ** argv) {
    EigenSolver<grid::RealScalarAntisymmetricPhysics> eigenSolver;
    double gridSpacing = 0.5/pow(2,2);
    double bounds = 20;
    Hamiltonian hamiltonian(gridSpacing, bounds);
    cout << "Grid spacing: " << gridSpacing << "\n";
    cout << "Grid points: " << hamiltonian.totalGridPoints << "\n";

    auto matvec = [&hamiltonian](vector<double> & ret, vector<double> const & v) {
        hamiltonian.matvec(ret, v);
        grid::RealScalarAntisymmetricPhysics::antiSymmetrize(ret, hamiltonian.numGridPointsPerDimension);
    };


    ofstream forceOfs("/tmp/forceMP");
    ofstream totalEnergyOfs("/tmp/totalEnergyMP");
    ofstream kineticEnergyOfs("/tmp/kenergyMP");
    ofstream ionElectronEnergyOfs("/tmp/ieenergyMP");
    ofstream electronElectrionEnergyOfs("/tmp/eeenergyMP");
    //while(ion1R < -0.25) {
        //ion1R += 0.05;
        //ion2R = -ion1R;
        if(ion1R == 0) {
            //break;
        }
        pair<double, double> spectrumBounds = eigenSolver.spectrumBounds(matvec, hamiltonian.numGridPointsPerDimension);
        cout << "Spectrum bounds: " << spectrumBounds.first << " " << spectrumBounds.second << "\n";

#if 0
        vector<vector<double>> subspace = eigenSolver.generateSubspace(matvec, std::min(hamiltonian.totalGridPoints, 500), hamiltonian.numGridPointsPerDimension);
#else
        vector<vector<double>> subspace = eigenSolver.generateTriangleFilteredSubspace(matvec,
                                                                                       std::min(hamiltonian.totalGridPoints, 100),
                                                                                       1,
                                                                                       hamiltonian.numGridPointsPerDimension,
                                                                                       spectrumBounds,
                                                                                       make_pair(spectrumBounds.first -
                                                                                                 5.0,
                                                                                                 spectrumBounds.first +
                                                                                                 5.0), 30);
#endif

        cout << "Starting Rayleigh-Ritz" << endl;
        EigenSolver<grid::RealScalarAntisymmetricPhysics>::EigenSystem eigenSystem = eigenSolver.rayleighRitz(matvec, subspace);
        cout << "extremal eigenvalues: " << eigenSystem.eigenvalues.front() << " " << eigenSystem.eigenvalues.back()
             << "\n";
        int numOutputFiles = 20;
        /*
         * Write 1D integrated density
         */
        {
            for(int k = 0; k < numOutputFiles; ++k) {
                stringstream sstr;
                sstr << "/tmp/densityMP" << k;
                ofstream ofs(sstr.str());
                vector<double> density(hamiltonian.numGridPointsPerDimension);
                for (int i = 0; i < hamiltonian.numGridPointsPerDimension; ++i) {
                    double v = 0;
                    for (int j = 0; j < hamiltonian.numGridPointsPerDimension; ++j) {
                        double t = grid::RealScalarAntisymmetricPhysics::elem(eigenSystem.eigenvectors[k], i, j, hamiltonian.numGridPointsPerDimension);
                        /*
                         * Density is normalized to integrate to 1, multiply by 2 to get the equivalent single particle density
                         */
                        v += 2 * t * t / gridSpacing;
                    }
                    ofs << hamiltonian.gridPoints[i] << " " << v << "\n";
                    density[i] = v;
                }
                stringstream sstr2;
                sstr2 << "/tmp/hartreePotentialMP" << k;
                ofstream ofs2(sstr2.str());
                vector<double> hartreePotential(hamiltonian.numGridPointsPerDimension);
                for (int i = 0; i < hamiltonian.numGridPointsPerDimension; ++i) {
                    double v = 0;
                    for (int j = 0; j < hamiltonian.numGridPointsPerDimension; ++j) {
                        v += -density[j] * hamiltonian.potential(hamiltonian.gridPoints[j], 1, hamiltonian.gridPoints[i]) * gridSpacing;
                    }
                    hartreePotential[i] = v;
                    ofs2 << hamiltonian.gridPoints[i] << " " << v << "\n";
                }
                stringstream sstr3;
                sstr3 << "/tmp/eePotentialMP" << k;
                ofstream ofs3(sstr3.str());
                stringstream sstr4;
                vector<double> effectiveSingleParticleEEPotential(hamiltonian.numGridPointsPerDimension);
                vector<double> vxc(hamiltonian.numGridPointsPerDimension);
                vector<pair<double,double>> densityVxc;
                for (int i = 0; i < hamiltonian.numGridPointsPerDimension; ++i) {
                    double v = 0;
                    for (int j = 0; j < hamiltonian.numGridPointsPerDimension; ++j) {
                        double t = grid::RealScalarAntisymmetricPhysics::elem(eigenSystem.eigenvectors[k], i, j, hamiltonian.numGridPointsPerDimension);
                        v += -t*t * hamiltonian.potential(hamiltonian.gridPoints[j], 1, hamiltonian.gridPoints[i]) / gridSpacing;
                    }
                    effectiveSingleParticleEEPotential[i] = v/density[i];
                    vxc[i] = v/density[i] - hartreePotential[i]/2;
                    ofs3 << hamiltonian.gridPoints[i] << " " << v/density[i] << "\n";
                    densityVxc.emplace_back(density[i], vxc[i]);
                }
                if(k == 0) {
                    ofstream ofs4("/tmp/vxcMP");
                    sort(begin(densityVxc), end(densityVxc), [](auto & x, auto & y) {
                        return x.first < y.first;
                    });
                    double lastDensity = -10;
                    for(auto & p : densityVxc) {
                        if(p.first > 1E-16) {
                            if(fabs(lastDensity - p.first)/p.first > 1E-3) {
                                ofs4 << p.first << " " << p.second << "\n";
                                lastDensity = p.first;
                            }
                        }
                    }
                }
                if(k == 0) {
                    double groundStateHartreeEnergy = 0;
                    double groundStateExcPlusHartreeEnergy = 0;
                    for(int i = 0; i < hamiltonian.numGridPointsPerDimension; ++i) {
                        groundStateHartreeEnergy += density[i] * hartreePotential[i] * gridSpacing;
                        groundStateExcPlusHartreeEnergy += density[i] * (vxc[i] + hartreePotential[i]) * gridSpacing;
                    }
                    cout << "Ground state single-particle equivalent Hartree energy: " << groundStateHartreeEnergy << "\n";
                    cout << "Ground state single-particle equivalent e-e energy: " << groundStateExcPlusHartreeEnergy << "\n";
                }
            }
        }
        /*
         * Write 2D orbitals
         */
        for (int k = 0; k < numOutputFiles; ++k) {
            stringstream sstr;
            sstr << "/tmp/eig" << k;
            double antiSymm = antiSymmetricMeasure(eigenSystem.eigenvectors[k], hamiltonian.numGridPointsPerDimension);
            double symm = symmetricMeasure(eigenSystem.eigenvectors[k], hamiltonian.numGridPointsPerDimension);
            string symmStr;
            if(antiSymm < 1E-5) {
                symmStr = "anti-symmetric";
            } else if(symm < 1E-5) {
                symmStr = "symmetric";
            } else {
                symmStr = "non-symmetric";
            }
            cout << k << " " << eigenSystem.eigenvalues[k] << " " << symmStr << "\n";
            ofstream ofs(sstr.str());
            for (int i = 0; i < hamiltonian.numGridPointsPerDimension; ++i) {
                for (int j = 0; j < hamiltonian.numGridPointsPerDimension; ++j) {
                    double t = grid::RealScalarAntisymmetricPhysics::elem(eigenSystem.eigenvectors[k], i, j, hamiltonian.numGridPointsPerDimension);
                    ofs << hamiltonian.gridPoints[i] << " " << hamiltonian.gridPoints[j] << " " << t << "\n";
                }
            }
        }
        /*
         * Compute "atomic" radii, for one atom
         */
        {
            int numAtoms = 2;
            vector<pair<double, double>> distanceDensity;
            for (int i = 0; i < hamiltonian.numGridPointsPerDimension; ++i) {
                for (int j = 0; j < hamiltonian.numGridPointsPerDimension; ++j) {
                    double x = hamiltonian.gridPoints[i] - ion1R;
                    double y = hamiltonian.gridPoints[j] - ion2R;
                    distanceDensity.emplace_back(sqrt(x * x + y * y),
                                                 grid::RealScalarAntisymmetricPhysics::elem(eigenSystem.eigenvectors[0], i, j,
                                                         hamiltonian.numGridPointsPerDimension));
                }
            }
            sort(begin(distanceDensity), end(distanceDensity), [](auto &x, auto &y) {
                return x.first < y.first;
            });
            double totalDensity = 0;
            for (int i = 0; i < distanceDensity.size(); ++i) {
                totalDensity += distanceDensity[i].second * distanceDensity[i].second;
                if (totalDensity >= 0.99/(numAtoms == 1 ? 2 : 1)) {
                    cout << "Radius: " << distanceDensity[i].first << "\n";
                    break;
                }
            }
        }
        /*
         * Compute forces
         */
        double ion1ElecForce = 0;
        double ion2ElecForce = 0;
        for (int i = 0; i < hamiltonian.numGridPointsPerDimension; ++i) {
            for (int j = 0; j < hamiltonian.numGridPointsPerDimension; ++j) {
                double t1 = grid::RealScalarAntisymmetricPhysics::elem(eigenSystem.eigenvectors[0], i, j, hamiltonian.numGridPointsPerDimension);
                ion1ElecForce += -hamiltonian.force(hamiltonian.gridPoints[i], 1, ion1R) * (t1 * t1);
                ion1ElecForce += -hamiltonian.force(hamiltonian.gridPoints[j], 1, ion1R) * (t1 * t1);
                ion2ElecForce += -hamiltonian.force(hamiltonian.gridPoints[i], 1, ion2R) * (t1 * t1);
                ion2ElecForce += -hamiltonian.force(hamiltonian.gridPoints[j], 1, ion2R) * (t1 * t1);
            }
        }
        //ion1ElecForce *= 2;
        //ion2ElecForce *= 2;
        cout << "ion force: " << hamiltonian.ionForce(ion2R, 1, ion1R) << "\n";
        cout << "ion force: " << hamiltonian.ionForce(ion1R, 1, ion2R) << "\n";
        cout << "elec force: " << ion1ElecForce << "\n";
        cout << "elec force: " << ion2ElecForce << "\n";
        cout << "total force 1: " << hamiltonian.ionForce(ion2R, 1, ion1R) + ion1ElecForce << "\n";
        cout << "total force 2: " << hamiltonian.ionForce(ion1R, 1, ion2R) + ion2ElecForce << "\n";
        cout << "ion1R: " << ion1R << "\n";
        forceOfs << ion1R << " " << hamiltonian.ionForce(ion2R, 1, ion1R) + ion1ElecForce << "\n";
        Energy energy = hamiltonian.getEnergies(matvec, eigenSystem);
        totalEnergyOfs << ion1R << " " << energy.totalEnergy << "\n";
        kineticEnergyOfs << ion1R << " " << energy.kineticEnergy << "\n";
        ionElectronEnergyOfs << ion1R << " " << energy.ionElectronEnergy<< "\n";
        electronElectrionEnergyOfs << ion1R << " " << energy.electronElectronEnergy<< "\n";
        cout << "Total energy: " << energy.totalEnergy << "\n";
        cout << "kinetic energy: " << energy.kineticEnergy << "\n";
        cout << "ion-elec energy: " << energy.ionElectronEnergy<< "\n";
        cout << "ion-ion energy: " << energy.ionIonEnergy<< "\n";
        cout << "elec-elec energy: " << energy.electronElectronEnergy << "\n";
    //}
}
