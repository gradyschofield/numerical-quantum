//
// Created by grady on 4/26/21.
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
#include<MathUtil.h>

using namespace std;

/*
 *  This program calculates the same thing as SecondIntegral_a but uses the spherical symmetry about q to remove
 *  two degrees of freedom (q_theta, q_phi).
 */

typedef BigSum AccumulatorType;

int main(int argc, char ** argv) {
    long N = CommandLine::parseNumPoints(argc, argv, 1E6);
    Bounds bounds2(0, 2, 0, 0, 0, 0);
    Bounds bounds1(0, 1, 0, M_PI, 0, 2 * M_PI);
    int numThreads = thread::hardware_concurrency();
    N /= numThreads;
    vector<PartialWork<AccumulatorType>> partialWorkArray(numThreads);
    vector<thread> threads;
    timespec start = Time::startTimer();
    for(int t = 0; t < numThreads; ++t) {
        threads.emplace_back([&,t]() {
            Generator generator1(bounds1);
            Generator generator2(bounds2);
            PartialWork<AccumulatorType> & partialWork = partialWorkArray[t];
            for (int i = 0; i < N; ++i) {
                SphericalPoint k = generator1.generate();
                SphericalPoint p = generator1.generate();
                SphericalPoint q = generator2.generate();
                double volumeProduct = k.volumeElementAt() * p.volumeElementAt() * square(q.getR());
                if (k.getR() > 1 || p.getR() > 1 ||
                    SphericalPoint::lengthOfSum(k, q) < 1 ||
                    SphericalPoint::lengthOfSum(p, q) < 1 ) {

                    partialWork.addOutVolume(volumeProduct);
                } else {
                    partialWork.addInVolume(volumeProduct);
                    double f1 = square(SphericalPoint::lengthOfSum({p, k, q}));
                    double f2 = square(q.getR()) + q.dot(k + p);
                    partialWork.addIntegral(volumeProduct / (q.getR() * q.getR() * f1 * f2));
                }
            }
        });
    }
    for(thread & t : threads) {
        t.join();
    }
    double runtime = Time::stopTimer(start);
    PartialWork<AccumulatorType> result = PartialWork<AccumulatorType>::accumulate(partialWorkArray);
    //double volumeNormalizer = square(bounds1.getVolume()) * bounds2.getVolume() / (result.getInVolume() + result.getOutVolume());
    double volumeNormalizer = square(bounds1.getVolume()) * 8/3 / (result.getInVolume() + result.getOutVolume());
    cout << "Number of threads used: " << numThreads << "\n";
    cout << "Run time: " << runtime/3600 << " hours \n";
    cout << "Number of points (millions): " << (numThreads*N)/1E6 << "\n";
    cout << "Integration region volume by monte carlo: " << 4 * M_PI * volumeNormalizer * result.getInVolume() << "\n";
    cout << "Monte Carlo integral: " << 4 * M_PI * volumeNormalizer * result.getIntegral() << "\n";
    cout << "Exact integral: ???"<< "\n";
    return 0;
}
