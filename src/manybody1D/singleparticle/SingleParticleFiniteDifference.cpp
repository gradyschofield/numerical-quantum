//
// Created by grady on 7/19/20.
//

#include<grid/RealScalarPhysics.h>
#include<EigenSolver.h>
#include<manybody1D/singleparticle/ForceCurve.h>
#include<manybody1D/singleparticle/Snapshot.h>
#include<manybody1D/singleparticle/Hamiltonian.h>

int main(int argc, char ** argv) {
    EigenSolver<grid::RealScalarPhysics> eigenSolver;
    double gridSpacing = 0.5/pow(2,2);
    double bounds = 20;
    double ion1R = -0.6;
    double ion2R = -ion1R;

    Hamiltonian hamiltonian(gridSpacing, bounds, ion1R, ion2R, "hartree");
    cout << "Grid spacing: " << gridSpacing << "\n";
    cout << "Grid points: " << hamiltonian.ndim << "\n";

    auto matvec = [&hamiltonian](vector<double> & ret, vector<double> const & v) {
        hamiltonian.matvec(ret, v);
    };

    /*
     * forceCurve
     * snapshot
     *
     * optional:
     * hartree
     * orbital dependent hartree
     * hartree plus correlation
     */
    string calculationType = "forceCurve";
    if(calculationType == "forceCurve") {
        ForceCurve::compute(true, ion1R, ion2R, hamiltonian, matvec, eigenSolver);
    } else if(calculationType == "snapshot") {
        Snapshot::compute(true, ion1R, ion2R, hamiltonian, matvec, eigenSolver);
    }

    if(hamiltonian.hasPrintVxcRangeWarning) {
        cout << "******** WARNING: The density has gone out of range of the Vxc functional *********\n";
    }
}