//
// Created by grady on 7/19/20.
//

#ifndef ATOM_CHEBYSHEVFILTER_H
#define ATOM_CHEBYSHEVFILTER_H

#include<vector>
#include<cmath>
#include<fstream>
#include<iostream>

using namespace std;


class ChebyshevFilter {
public:
    vector<double> coef;
    pair<double, double> spectrumBounds;

    static vector<double> chebyshevNodes(int n) {
        vector<double>ret(n);
        for(int i = 1; i <= n; ++i) {
            ret[i-1] = cos(M_PI*(2*i-1)/(2*n));
        }
        return ret;
    }

    static double evaluateChebyshev(int p, double x) {
        if(p == 0) {
            return 1/sqrt(M_PI);
        } else if(p == 1) {
            return x/sqrt(M_PI_2);
        } else {
            double c0 = 1;
            double c1 = x;
            double c2;
            for (int i = 2; i <= p; ++i) {
                c2 = 2 * x * c1 - c0;
                c0 = c1;
                c1 = c2;
            }
            return c2/sqrt(M_PI_2);
        }
    }

    static vector<double> evaluateChebyshev(int p, vector<double> const & nodes) {
        vector<double> ret;
        ret.reserve(nodes.size());
        for(double x : nodes) {
            ret.push_back(evaluateChebyshev(p, x));
        }
        return ret;
    }

    static ChebyshevFilter triangle(int n, double low, double high, double start, double end) {
        auto filter = [low,high,start,end](double x) {
            double shiftedStart = 2 * (start - low) / (high - low) - 1;
            double shiftedEnd = 2 * (end - low) / (high - low) - 1;
            double mid = (shiftedStart + shiftedEnd)/2;
            if(x < shiftedStart) {
                return 0.0;
            } else if(x < mid) {
                return (x - shiftedStart) / (mid-shiftedStart);
            } else if(x< shiftedEnd) {
                return (shiftedEnd - x) / (shiftedEnd-mid);
            } else {
                return 0.0;
            }
        };
        return fromFilter(filter, n, low, high, true);
    }

    static ChebyshevFilter step(int n, double low, double high, double start, double end) {
        auto filter = [&](double x) {
            double shiftedStart = 2 * (start - low) / (high - low) - 1;
            double shiftedEnd = 2 * (end - low) / (high - low) - 1;
            if(x< shiftedStart) {
                return 0.0;
            } else if(x < shiftedEnd) {
                return 1.0;
            } else {
                return 0.0;
            }
        };
        return fromFilter(filter, n,  low, high, true);
    }

    template<typename Filter>
    static ChebyshevFilter fromFilter(Filter && filter, int n, double low, double high, bool smoothed = false) {
        vector<double> nodes = chebyshevNodes(n);
        vector<double> coef;
        for(int p = 0; p < n; ++p) {
            ofstream ofs("/tmp/filt");

            vector<double> chebp = evaluateChebyshev(p, nodes);
            if(n == 0) {

            } else {
                double d = 0;
                for(int i = 0; i < nodes.size(); ++i) {
                    ofs << nodes[i] << " " << filter(nodes[i]) << "\n";
                    if(p > 0) {
                        d += M_PI / n * filter(nodes[i]) * chebp[i];
                    } else {
                        d += M_PI / n * filter(nodes[i]) * chebp[i];
                    }
                }
                coef.push_back(d);
            }
        }
        if(smoothed) {
            double k = coef.size();
            double alpha = M_PI / (k + 2);
            for (int i = 0; i < coef.size(); ++i) {
                coef[i] *= ((1 - i / (k + 2)) * sin(alpha) * cos(i * alpha) + 1 / (k + 2) * cos(alpha) * sin(i * alpha)) / sin(alpha);
            }
        }
        return ChebyshevFilter{move(coef), make_pair(low, high)};
    }

    void plot(string filename) {
        int degree = coef.size()-1;
        ofstream ofs(filename);
        for(double x = -1; x < 0; x += 0.001) {
            double y = 0;
            for(int p = 0; p <= degree; ++p) {
                y += coef[p] * evaluateChebyshev(p, x);
            }
            ofs << x << " " << y << "\n";
        }
        double y = 0;
        for(int p = 0; p <= degree; ++p) {
            y += coef[p] * evaluateChebyshev(p, 0);
        }
        ofs << 0 << " " << y << "\n";
        for(double x = 0.001; x < 1; x += 0.001) {
            double y = 0;
            for(int p = 0; p <= degree; ++p) {
                y += coef[p] * evaluateChebyshev(p, x);
            }
            ofs << x << " " << y << "\n";
        }
    }
};


#endif //ATOM_CHEBYSHEVFILTER_H
