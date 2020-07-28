//
// Created by grady on 7/23/20.
//

#ifndef ATOM_SNAPSHOT_H
#define ATOM_SNAPSHOT_H

#include<fstream>
#include<iostream>
#include<vector>
#include<manybody1D/singleparticle/Hamiltonian.h>
#include<EigenSolver.h>
#include<grid/RealScalarPhysics.h>

using namespace std;

class Snapshot {
public:
    template<typename Matvec>
    static void compute(bool chebyshevFilter, double ion1R, double ion2R, Hamiltonian & hamiltonian,
                        Matvec &&matvec, EigenSolver<grid::RealScalarPhysics> const &eigenSolver) {
        hamiltonian.setIon1R(ion1R);
        hamiltonian.setIon2R(ion2R);
        pair<double, double> spectrumBounds = eigenSolver.spectrumBounds(matvec, hamiltonian.ndim);
        cout << "Spectrum bounds: " << spectrumBounds.first << " " << spectrumBounds.second << "\n";

        hamiltonian.startCalculation();
        EigenSolver<grid::RealScalarPhysics>::EigenSystem eigenSystem;
        do {
            vector<vector<double>> subspace;
            if (!chebyshevFilter) {
                subspace = eigenSolver.generateSubspace(matvec, std::min(hamiltonian.ndim, 500),
                                                        hamiltonian.ndim);
            } else {
                subspace = eigenSolver.generateTriangleFilteredSubspace(matvec,
                                                                        std::min(hamiltonian.ndim, 100),
                                                                        1,
                                                                        hamiltonian.ndim,
                                                                        spectrumBounds,
                                                                        make_pair(spectrumBounds.first - 2.5,
                                                                                  spectrumBounds.first + 2.5),
                                                                        30);
            }

            cout << "Starting Rayleigh-Ritz" << endl;
            eigenSystem = eigenSolver.rayleighRitz(matvec, subspace);
            hamiltonian.updatePotential(eigenSystem);
        } while(!hamiltonian.isConverged());
        cout << "extremal eigenvalues: " << eigenSystem.eigenvalues.front() << " " << eigenSystem.eigenvalues.back()
             << "\n";
        for (int j = 0; j < 10; ++j) {
            stringstream sstr;
            sstr << "/tmp/eig" << j;
            cout << j << " " << eigenSystem.eigenvalues[j] << "\n";
            ofstream ofs(sstr.str());
            for (int i = 0; i < hamiltonian.gridPoints.size(); ++i) {
                ofs << hamiltonian.gridPoints[i] << " " << eigenSystem.eigenvectors[j][i] << "\n";
            }
            if(j==0) {
                ofstream ofs("/tmp/groundStateDensitySP");
                for (int i = 0; i < hamiltonian.gridPoints.size(); ++i) {
                    double t = eigenSystem.eigenvectors[j][i];
                    ofs << hamiltonian.gridPoints[i] << " " << 2*t*t/hamiltonian.h<< "\n";
                }
            }
        }
        double ion1ElecForce = 0;
        double ion2ElecForce = 0;
        for (int i = 0; i < hamiltonian.gridPoints.size(); ++i) {
            double t1 = eigenSystem.eigenvectors[0][i];
            double t2 = eigenSystem.eigenvectors[0][i];
            ion1ElecForce += -hamiltonian.force(hamiltonian.gridPoints[i], 1, ion1R) * (t1 * t1 + t2 * t2);
            ion2ElecForce += -hamiltonian.force(hamiltonian.gridPoints[i], 1, ion2R) * (t1 * t1 + t2 * t2);
        }
        cout << "ion force: " << hamiltonian.ionForce(ion2R, 1, ion1R) << "\n";
        cout << "ion force: " << hamiltonian.ionForce(ion1R, 1, ion2R) << "\n";
        cout << "elec force: " << ion1ElecForce << "\n";
        cout << "elec force: " << ion2ElecForce << "\n";
        cout << "total force 1: " << hamiltonian.ionForce(ion2R, 1, ion1R) + ion1ElecForce << "\n";
        cout << "total force 2: " << hamiltonian.ionForce(ion1R, 1, ion2R) + ion2ElecForce << "\n";
        Energy energy = hamiltonian.getEnergies(matvec, eigenSystem);
        cout << "Total energy: " << energy.totalEnergy << "\n";
        cout << "kinetic energy: " << energy.kineticEnergy << "\n";
        cout << "ion-elec energy: " << energy.ionElectronEnergy<< "\n";
        cout << "ion-ion energy: " << energy.ionIonEnergy<< "\n";
        cout << "elec-elec energy: " << energy.electronElectronEnergy << "\n";
        cout << "xc energy: " << energy.exchangeCorrelationEnergy << "\n";
    }
};
#endif //ATOM_SNAPSHOT_H
