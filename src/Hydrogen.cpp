#include <iostream>
#include<cmath>
#include<lapacke.h>
#include<fstream>
#include<vector>

using namespace std;

/*
 * This program will compute radial wavefunctions for the hydrogen atom for
 * a given l quantum number using linear finite elements.  The radial quantum number, l,
 * boundary of the radial grid, and number eigenvectors to compute are defined below. The
 * function to get a variable grid spacing is good up to 10,000 a.u. where about 50
 * l=0 radial states are converged.  If you want to go further, you need to add new cases
 * to the function that computes grid spacing.  Keep in mind that increasing the l quantum
 * number makes the states for a given radial quantum number more diffuse.
 *
 * The eigenvalues are print to the console and the eigenfunctions are written to files /tmp/eig*
 * The PlotEigs.py script will print the wavefunctions or the logarithm of the wavefunctions.
 *
 * It is good to look at the log of the eigenfunction to get a clear picture of how much space there is
 * between the boundary and where the last oscillations occur.  There needs to be plenty of space for
 * the eigenfunction to go to zero.  If there is not enough space, then what you will find if you make
 * the domain bigger is that new eigenvalues appear near the top of the discrete part of the spectrum where
 * you might have thought you had a converged spectrum.
 */

double l = 0;
int numEigenvectors = 100;
double boundaryEnd = 10000;

double getGridSpacing(double x) {
#if 1
    //First grid spacing is 0.05 au
    if(x < 400) {
        return 0.05 * exp(0.01 * x);
    } else if(x < 2000){
        return 2.5 * exp(0.001 *x);
    } else {
        return 18.5 * exp(0.000125 * x); // good up to about 10,000 au
    }
#else
    if(x < 400) {
        return h * exp(0.008 * x);
    } else {
        return 1.5 * exp(0.001 *x);
    }
#endif
}

void getEigenSystem(int numEigenvectors, vector<double> & stiffness, vector<double> & mass, vector<double> & eigenvalues, vector<double> & eigenvectors) {
    int itype = 1;
    char jobz = 'V';
    char uplo = 'U';
    char range = 'I';
    int n = (int)sqrt(stiffness.size());
    cout << n << "\n";
    int lda = n;
    int ldb = n;
    int ldz = n;
    int numEigenvaluesFound;
    eigenvalues.resize(n);
    eigenvectors.resize(n * n);
    vector<int>ifail(n);
    LAPACKE_dsygvx(LAPACK_COL_MAJOR, itype, jobz, range, uplo, n, stiffness.data(), lda, mass.data(), ldb, 0.0, 0.0,
            1, numEigenvectors, 0, &numEigenvaluesFound, eigenvalues.data(), eigenvectors.data(), ldz, ifail.data());
    cout << "Num eigenvalue found: " << numEigenvaluesFound << "\n";
}

int main() {
    vector<double> knots;
    double x = 0;
    vector<double> gridSpacing;
    while(x < boundaryEnd) {
        knots.push_back(x);
        double h = getGridSpacing(x);
        gridSpacing.push_back(h);
        x += h;
    }
    size_t numKnots = knots.size();
    cout << "Num knots: " << numKnots << "\n";
    cout << "Max grid spacing: " << gridSpacing.back() << "\n";
    gridSpacing.push_back(boundaryEnd - knots.back());
    vector<double> stiffness(numKnots*numKnots);
    vector<double> mass(numKnots*numKnots);
    int n = numKnots;
    for(int i = 0; i < numKnots; ++i) {
        double h1 = 0;
        if (i > 0) {
            h1 = gridSpacing[i - 1];
        }
        double h2 = gridSpacing[i];
        double r = knots[i];
        double diagonalGrad;
        double a = r - h1;
        double b = r;
        double c = r + h2;

        /*
         * Gradient contribution
         */
        double i1 = (pow(b, 3) - pow(a, 3))/(3*h1*h1);
        double i2 = (pow(c, 3) - pow(b, 3))/(3*h2*h2);
        if (i == 0) {
            diagonalGrad = 0.5*i2; //0.5 coming from Schrodinger
        } else {
            //diagonalGrad = 0.5*r*r*(1/h1+1/h2); //0.5 coming from Schrodinger
            diagonalGrad = 0.5*(i1+i2); //0.5 coming from Schrodinger
        }
        double offDiagonalGrad = -0.5*i2; //0.5 coming from Schrodinger


        /*
         * Angular momentum contributions
         */
        double diagonalAngular;
        if (i == 0) {
            diagonalAngular = 0.5*l*(l+1)*(h2/3);
        } else {
            diagonalAngular = 0.5*l*(l+1)*(h1/3+h2/3);
        }
        double offDiagonalAngular = 0.5*l*(l+1)*(h2/6);

        /*
         * Coulomb contribution
         */
        double diagonalCoulomb;
        auto coulombDiagonalIntegralExpr = [h2](double x, double off, double h) {
            return -(pow(x, 4)/4 - 2*off*pow(x,3)/3 + off*off*pow(x,2)/2) / (h*h);
        };
        i1 = coulombDiagonalIntegralExpr(b, a, h1) - coulombDiagonalIntegralExpr(a, a, h1);
        i2 = coulombDiagonalIntegralExpr(c, c, h2) - coulombDiagonalIntegralExpr(b, c, h2);
        if(i == 0) {
            diagonalCoulomb = i2;
        } else {
            diagonalCoulomb = i1+i2;
        }
        auto coulombOffDiagonalIntegralExpr = [h2](double x, double b, double c) {
            return -((c+b)/3*pow(x,3) - b*c/2*pow(x,2) - pow(x,4)/4)/(h2*h2);
        };
        i2 = coulombOffDiagonalIntegralExpr(c, b, c) - coulombOffDiagonalIntegralExpr(b, b, c);
        double offDiagonalCoulomb = i2;

        /*
         * Mass matrix elements
         */
        auto massDiagonalExpr = [](double x, double off, double h) {
            return (pow(x, 5)/5 - off*pow(x,4)/2 + off*off*pow(x,3)/3) / (h*h);
        };
        i1 = massDiagonalExpr(b, a, h1) - massDiagonalExpr(a, a, h1);
        i2 = massDiagonalExpr(c, c, h2) - massDiagonalExpr(b, c, h2);
        double diagonalMass;
        if(i == 0) {
            diagonalMass = i2;
        } else {
            diagonalMass = i1 + i2;
        }
        auto massOffDiagonalExpr = [h2](double x, double b, double c) {
            return ((c+b)/4*pow(x,4) - b*c/3*pow(x,3) - pow(x,5)/5)/(h2*h2);
        };
        i2 = massOffDiagonalExpr(c, b, c) - massOffDiagonalExpr(b, b, c);
        double offDiagonalMass = i2;

        stiffness[i+n*i] += diagonalGrad + diagonalAngular + diagonalCoulomb;
        if(i < numKnots-1) {
            stiffness[i + n * (i + 1)] += offDiagonalGrad + offDiagonalAngular + offDiagonalCoulomb;
        }

        mass[i+n*i] += diagonalMass;
        if(i < numKnots-1) {
            mass[i + n * (i + 1)] += offDiagonalMass;
        }
    }
    vector<double> eigenvalues;
    vector<double> eigenvectors;
    cout << "Compute system" << std::endl;
    getEigenSystem(numEigenvectors, stiffness, mass, eigenvalues, eigenvectors);
    for(int i = 0; i < numEigenvectors; ++i) {
        cout << (i+1) << " " << eigenvalues[i] << "\n";
    }
    for(int j = 0; j < numEigenvectors; ++j) {
        stringstream sstr;
        sstr << "/tmp/eig" << (j+1);
        ofstream ofs(sstr.str());
        for (int i = 0; i < numKnots; ++i) {
            ofs << knots[i] << " " << eigenvectors[i + j * numKnots] << "\n";
        }
    }

    return 0;
}
