//
// Created by grady on 7/23/20.
//

#ifndef ATOM_HAMILTONIAN_H
#define ATOM_HAMILTONIAN_H

#include<Constants.h>
#include<stdeig/EigenSolver.h>
#include <Energy.h>
#include <algorithm>

class Hamiltonian {
public:
    vector<double> deriv2_2{1, -2, 1};
    vector<double> deriv2_4{-1.0/12, 4.0/3, -5.0/2, 4.0/3, -1.0/12};

    int ndim;
    double h;
    vector<double> deriv2coef;
    int deriv2Length;
    int deriv2BackOffset;
    vector<double> gridPoints;
    double ion1R;
    double ion2R;
    bool needsSCF;
    string potentialType;
    vector<double> hartreePotential;
    vector<vector<double>> orbitalDependentHartreePotential;
    double scfResidual = std::numeric_limits<double>::max();
    vector<pair<double, double>> vxc;
    vector<pair<double, double>> exc;
    vector<double> density;
    double mixingParam = 0.3;

    void setIon1R(double x) {
        ion1R = x;
    }

    void setIon2R(double x) {
        ion2R = x;
    }

    void startCalculation() {
        ofstream ofs("/tmp/ionpot");
        for (double x : gridPoints) {
            ofs << x << " " << potential(x, 1, ion1R) + potential(x, 1, ion2R) << "\n";
        }
        if(!needsSCF) {
            return;
        }
        /*
         * Reset potentials to zero?
         */
    }

    bool hasPrintVxcRangeWarning = false;

    double interpolateXc(double d, vector<pair<double,double>> const & func) {
        if(d < func[0].first) {
            return func[0].second;
        }
        for(int i = 0; i < func.size(); ++i) {
            if(func[i].first > d) {
                double t = (d- func[i-1].first) / (func[i].first - func[i-1].first);
                return t * func[i].second + (1-t) * func[i-1].second;
            }
        }
        if(!hasPrintVxcRangeWarning) {
            cout << "******** WARNING: The density has gone out of range of the Vxc functional *********\n";
            hasPrintVxcRangeWarning = true;
        }
        return func.back().second;
    }

    Hamiltonian(double gridSpacing, double bounds, double ion1R, double ion2R, string potentialType) {
        h = gridSpacing;
        ndim = 2 * bounds / gridSpacing;
        for (int i = ndim / 2 - 1; i >= 0; --i) {
            gridPoints.push_back(-h / 2 - i * h);
        }
        for (int i = 0; i < ndim / 2; ++i) {
            gridPoints.push_back(h / 2 + i * h);
        }
        ndim = gridPoints.size();
        deriv2coef = deriv2_2;
        for (auto &coef : deriv2coef) {
            coef *= -1 / (2 * constants::massOfElectron * h * h);
        }
        for (auto &coef : deriv2coef) {
            cout << "coef: " << coef << "\n";
        }
        deriv2Length = deriv2coef.size();
        deriv2BackOffset = (deriv2Length - 1) / 2;
        ofstream ofs("/tmp/ionpot");
        for (double x : gridPoints) {
            ofs << x << " " << potential(x, 1, ion1R) + potential(x, 1, ion2R) << "\n";
        }
        Hamiltonian::ion1R = ion1R;
        Hamiltonian::ion2R = ion2R;
        Hamiltonian::potentialType = potentialType;
        if (potentialType == "noInteraction") {
            needsSCF = false;
        } else if (potentialType == "hartree") {
            needsSCF = true;
            density = vector<double>(ndim);
            hartreePotential = vector<double>(ndim);
        } else if (potentialType == "orbitalProjectedHartree") {
            needsSCF = true;
            for(int i = 0; i < 2; ++i) {
                orbitalDependentHartreePotential[i] = vector<double>(ndim);
            }
        }
        ifstream ifs("/tmp/vxcMP");
        ofstream computedExcOfs("/tmp/computedExc");
        while(true) {
            string line;
            getline(ifs, line);
            if(ifs.fail()) {
                break;
            }
            stringstream  str(line);
            double d, potential;
            str >> d >> potential;
            exc.emplace_back(d, d * potential);
            computedExcOfs << d << " " << d*potential << "\n";
        }
        ofstream computedVxcOfs("/tmp/computedVxc");
        double v = 0;
        for(int i = 1; i < exc.size(); ++i) {
            //v += 0.5 * (vxc[i].second + vxc[i-1].second)*(vxc[i].first - vxc[i-1].first);
            //double x = (vxc[i].first + vxc[i-1].first)*0.5;
            double v = (exc[i].second - exc[i-1].second)/(exc[i].first - exc[i-1].first);
            double x = exc[i].first;
            vxc.emplace_back(x, v);
            computedVxcOfs << x << " " << v << "\n";
        }
    }

    void updatePotential(EigenSolver<grid::RealScalarPhysics>::EigenSystem const & eigenSystem) {
        density = vector<double>(ndim);
        for(int i = 0; i < 2; ++i) {
            for(int j = 0; j < ndim; ++j) {
                double t = eigenSystem.eigenvectors[0][j];
                density[j] += t*t / h;
            }
        }
        vector<double> pot = electronPotential(density);
        ofstream vofs("/tmp/epot");
        ofstream dofs("/tmp/density");
        double n = 0;
        for(int i = 0; i < ndim; ++i) {
            vofs << gridPoints[i] << " " << pot[i] << "\n";
            dofs << gridPoints[i] << " " << density[i] << "\n";
            n += density[i] * h;
        }
        cout << "density integral: " << n << "\n";
        if(!needsSCF) {
            return;
        }
        scfResidual = 0;
        if(potentialType == "hartree") {
            for(int i = 0; i < ndim; ++i) {
                scfResidual += fabs(hartreePotential[i] - pot[i]) * density[i] * h;
                hartreePotential[i] = (1-mixingParam) * hartreePotential[i] + mixingParam * pot[i];
            }
        }
        cout << "scf residual: " << scfResidual << "\n";
    }

    bool isConverged() {
        if(!needsSCF) {
            return true;
        }
        return scfResidual < 1E-5;
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
     * v is a wfn squared or the sum of wfns squared
     */
    vector<double> electronPotential(vector<double> const & v) {
        vector<double> pot(ndim);
        for(int i = 0; i < ndim; ++i) {
            for(int j = 0; j < ndim; ++j) {
                pot[i] += -v[j]*potential(gridPoints[j], 1, gridPoints[i]) * h;
            }
        }
        return pot;
    }

    /*
     *  This is the force for like charges.  Use -force when integrating election-ion contributions.
     */
    double force(double x, double N, double R) const {
        double t = fabs(x - R);
        double sign = (x-R) < 0 ? -1 : 1;
        if(t < 0.5) {
            return -sign*constants::chargeOfElectron*N*(8*t);
        } else {
            return -sign*constants::chargeOfElectron*N/(t*t);
        }
    }

    double ionForce(double x, double N, double R) const {
        double t = fabs(x - R);
        double sign = (x-R) < 0 ? -1 : 1;
        return -sign*constants::chargeOfElectron*N/(t*t);
    }

    void matvecKinetic(vector<double> & ret, vector<double> const & v) {
        for(int i = 0; i < v.size(); ++i) {
            for(int j = std::max(0, i-deriv2BackOffset); j <= std::min((int)v.size()-1, i + deriv2BackOffset); ++j) {
                ret[i] += deriv2coef[j - i + deriv2BackOffset] * v[j];
            }
        }
    }

    void matvecIonElectron(vector<double> & ret, vector<double> const & v) {
        for(int i = 0; i < v.size(); ++i) {
            ret[i] += v[i] * (potential(gridPoints[i], 1, ion2R) + potential(gridPoints[i], 1, ion1R));
        }
    }

    void matvecElectronElectron (vector<double> & ret, vector<double> const & v) {
        if(potentialType == "hartree") {
            for (int i = 0; i < ndim; ++i) {
                ret[i] += v[i] * hartreePotential[i];
            }
        }
    }

    void matvecExchangeCorrelaction(vector<double> & ret, vector<double> const & v) {
        if(potentialType == "hartree") {
            for (int i = 0; i < ndim; ++i) {
                ret[i] += v[i] * interpolateXc(density[i], vxc);
            }
        }
    }

    double energyExchangeCorrelation() {
        double e = 0;
        if(potentialType == "hartree") {
            for (int i = 0; i < ndim; ++i) {
                e += interpolateXc(density[i], exc) * h;
            }
        }
        return e;
    }

    void matvec(vector<double> & ret, vector<double> const & v) {
        for(auto & x : ret) x = 0;
        matvecKinetic(ret, v);
        matvecIonElectron(ret, v);
        matvecElectronElectron(ret, v);
        matvecExchangeCorrelaction(ret, v);
    }

    template<typename Matvec>
    Energy getEnergies(Matvec && matvec, EigenSolver<grid::RealScalarPhysics>::EigenSystem const & eigenSystem) {
        vector<double> ret(ndim);
        matvecKinetic(ret, eigenSystem.eigenvectors[0]);
        double kineticEnergy = 2 * grid::RealScalarPhysics::dot(eigenSystem.eigenvectors[0], ret);

        for(auto & x : ret) x = 0;
        matvecIonElectron(ret, eigenSystem.eigenvectors[0]);
        double ionElectronEnergy = 2 * grid::RealScalarPhysics::dot(eigenSystem.eigenvectors[0], ret);

        for(auto & x : ret) x = 0;
        matvecElectronElectron(ret, eigenSystem.eigenvectors[0]);
        double electronElectronEnergy = grid::RealScalarPhysics::dot(eigenSystem.eigenvectors[0], ret);

        double exchangeCorrelationEnergy = energyExchangeCorrelation();

        double ionIonEnergy = -potential(ion1R, 1, ion2R);
        return Energy(kineticEnergy, ionElectronEnergy, electronElectronEnergy, exchangeCorrelationEnergy, ionIonEnergy);
    }
};

#endif //ATOM_HAMILTONIAN_H
