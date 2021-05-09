//
// Created by Grady Schofield on 4/23/21.
//

#include<iostream>
#include<vector>
#include<thread>

#include<BigSum.h>
#include<Bounds.h>
#include<CommandLine.h>
#include<Generator.h>
#include<PartialWork.h>
#include<SphericalPoint.h>

using namespace std;

/*
 * This program calculates Integral( dk dq 1/q^2 Heaviside(1 - |k+q|) Heaviside(1 - |k|))
 * where k and q are 3-vectors.  This is part of the calculation of the first order
 * term in the perturbation expansion of the electron-electron energy for a homogenous
 * electron gas.  The full expression for this term is - e^2 4 pi V / (2 pi)^6 times this
 * integral.  Using V = 4 pi (r_0)^3 N / 3 and r_s = r_0 / a_0 (bohr radius), this can
 * be reduced to E(1)/N = e^2 / (2 a_0) (-0.916 / r_s).  The integral should be equal
 * to 4 pi^2.
 *
 * This integral can be done fairly easily by examining the geometry of the integrand
 * (see Fetter & Walecka p. 28), but integrals for higher terms in the perturbation
 * series must be done numerically. So this is a good place to get started with Monte
 * Carlo integration.
 */

typedef BigSum AccumulatorType;

int main(int argc, char ** argv) {
    long N = CommandLine::parseNumPoints(argc, argv, 1E6);
    Bounds bounds(0, 2, 0, M_PI, 0, 2 * M_PI);
    int numThreads = thread::hardware_concurrency();
    N /= numThreads;
    vector<PartialWork<AccumulatorType>> partialWorkArray(numThreads);
    vector<thread> threads;
    timespec t1, t2;
    clock_gettime(CLOCK_REALTIME, &t1);
    for(int t = 0; t < numThreads; ++t) {
        threads.emplace_back([&,t]() {
            Generator generator(bounds);
            PartialWork<AccumulatorType> & partialWork = partialWorkArray[t];
            for (int i = 0; i < N; ++i) {
                SphericalPoint k = generator.generate();
                SphericalPoint q = generator.generate();
                double volumeProduct = k.volumeElementAt() * q.volumeElementAt();
                if (k.getR() > 1 || SphericalPoint::lengthOfSum(k, q) > 1) {
                    partialWork.addOutVolume(volumeProduct);
                } else {
                    partialWork.addInVolume(volumeProduct);
                    partialWork.addIntegral(volumeProduct / (q.getR() * q.getR()));
                }
            }
        });
    }
    for(thread & t : threads) {
        t.join();
    }
    clock_gettime(CLOCK_REALTIME, &t2);
    double time = (t2.tv_sec * 1E9 + t2.tv_nsec - t1.tv_sec * 1E9 - t1.tv_nsec)/1E9;
    PartialWork<AccumulatorType> result = PartialWork<AccumulatorType>::accumulate(partialWorkArray);
    double volumeNormalizer = pow(bounds.getVolume(), 2) / (result.getInVolume() + result.getOutVolume());
    cout << "Number of threads used: " << numThreads << "\n";
    cout << "Number of points (millions): " << (numThreads*N)/1E6 << "\n";
    cout << "Integration region volume by monte carlo: " << volumeNormalizer * result.getInVolume() << "\n";
    double integral = volumeNormalizer * result.getIntegral();
    double exact =  4* M_PI * M_PI;
    cout << "Monte Carlo integral: " << integral << "\n";
    cout << "Exact integral: " << exact << "\n";
    cout << "Relative error: " << fabs(integral - exact)/exact << "\n";
    cout << "Time per sample " << time/(numThreads*N) << "\n";
    return 0;
}