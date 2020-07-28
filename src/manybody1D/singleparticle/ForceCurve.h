//
// Created by grady on 7/23/20.
//

#ifndef ATOM_FORCECURVE_H
#define ATOM_FORCECURVE_H

#include<fstream>
#include<iostream>
#include<vector>
#include<manybody1D/singleparticle/Hamiltonian.h>
#include<EigenSolver.h>
#include<grid/RealScalarPhysics.h>
#include <Energy.h>

using namespace std;

class ForceCurve {
public:
    template<typename Matvec>
    static void compute(bool chebyshevFilter, double ion1R, double ion2R, Hamiltonian & hamiltonian,
                        Matvec &&matvec, EigenSolver<grid::RealScalarPhysics> const &eigenSolver) {
        if (ion1R >= ion2R) {
            cout << "ForceCurve expects ion 1 to be on the left of ion 2\n";
            return;
        }
        ofstream forceOfs("/tmp/forceSP");
        ofstream totalEnergyOfs("/tmp/totalEnergySP");
        ofstream kineticEnergyOfs("/tmp/kenergySP");
        ofstream ionElectronEnergyOfs("/tmp/ieenergySP");
        ofstream electronElectrionEnergyOfs("/tmp/eeenergySP");
        while (true) {
            ion1R += 0.125/8.0;
            ion2R = -ion1R;
            //if (ion1R >= ion2R) {
            if (ion1R >= -0.2) {
                break;
            }
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
            cout << "extremal eigenvalues: " << eigenSystem.eigenvalues.front() << " " << eigenSystem.eigenvalues.back() << "\n";
            for (int j = 0; j < 10; ++j) {
                stringstream sstr;
                sstr << "/tmp/eig" << j;
                cout << j << " " << eigenSystem.eigenvalues[j] << "\n";
                ofstream ofs(sstr.str());
                for (int i = 0; i < hamiltonian.gridPoints.size(); ++i) {
                    ofs << hamiltonian.gridPoints[i] << " " << eigenSystem.eigenvectors[j][i] << "\n";
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
            cout << "force: " << hamiltonian.ionForce(ion2R, 1, ion1R) << "\n";
            cout << "force: " << hamiltonian.ionForce(ion1R, 1, ion2R) << "\n";
            cout << "e force: " << ion1ElecForce << "\n";
            cout << "e force: " << ion2ElecForce << "\n";
            cout << "total force 1: " << hamiltonian.ionForce(ion2R, 1, ion1R) + ion1ElecForce << "\n";
            cout << "total force 2: " << hamiltonian.ionForce(ion1R, 1, ion2R) + ion2ElecForce << "\n";
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
            cout << "xc energy: " << energy.exchangeCorrelationEnergy << "\n";
        }
    }
};


#endif //ATOM_FORCECURVE_H
