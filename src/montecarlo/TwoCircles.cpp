//
// Created by Grady Schofield on 4/23/21.
//
#include<random>
#include<thread>
#include<iostream>

using namespace std;

class Point {
    double x, y;
public:
    Point(double x, double y) : x(x), y(y) {}

    double getX() const {
        return x;
    }

    double getY() const {
        return y;
    }
};

class Circle {
    double x, y, r, r2;
public:
    Circle(double x, double y, double r)
        : x(x), y(y), r(r), r2(r*r)
    {
    }

    bool inside(Point const & p) {
        double dx = p.getX() - Circle::x;
        double dy = p.getY() - Circle::y;
        return dx*dx + dy*dy < r2;
    }

    static bool insideAll(Point const & p, vector<Circle> const & circles) {
        for(Circle const & c : circles) {
            double dx = p.getX() - c.getX();
            double dy = p.getY() - c.getY();
            if(dx*dx + dy*dy > c.getR2()) {
                return false;
            }
        }
        return true;
    }

    double getX() const {
        return x;
    }

    double getY() const {
        return y;
    }

    double getR() const {
        return r;
    }

    double getR2() const {
        return r2;
    }
};

class alignas(128) Bounds {
    double minx = numeric_limits<double>::max();
    double maxx = -numeric_limits<double>::max();
    double miny = numeric_limits<double>::max();
    double maxy = -numeric_limits<double>::max();

public:
    Bounds(){
    }

    Bounds(double minx, double maxx, double miny, double maxy)
        : minx(minx), maxx(maxx), miny(miny), maxy(maxy)
    {
    }

    void observed(Point const & p) {
        minx = min(minx, p.getX());
        maxx = max(maxx, p.getX());
        miny = min(miny, p.getY());
        maxy = max(maxy, p.getY());
    }

    static Bounds find(vector<Circle> const & circles) {
        double miny, minx = miny = numeric_limits<double>::max();
        double maxy, maxx = maxy = -numeric_limits<double>::max();
        for(Circle const & c : circles) {
            minx = min(minx, c.getX() - c.getR());
            maxx = max(maxx, c.getX() + c.getR());
            miny = min(miny, c.getY() - c.getR());
            maxy = max(maxy, c.getY() + c.getR());
        }
        return Bounds(minx, maxx, miny, maxy);
    }

    double getMinX() const {
        return minx;
    }

    double getMaxX() const {
        return maxx;
    }

    double getMinY() const {
        return miny;
    }

    double getMaxY() const {
        return maxy;
    }

    double getVolume() const {
        return (maxx - minx) * (maxy - miny);
    }

    Bounds grow(double p) {
        double midx = 0.5*(maxx + minx);
        double midy = 0.5*(maxy + miny);
        double widthx = 0.5*(maxx - minx);
        double widthy = 0.5*(maxy - miny);
        return Bounds(midx - widthx * p,
                      midx + widthx * p,
                      midy - widthy * p,
                      midy + widthy * p);
    }
};

class Generator {
    mt19937_64 generator;
    uniform_real_distribution<> unifx;
    uniform_real_distribution<> unify;
public:
    Generator(Bounds const & bounds) {
        random_device r;
        seed_seq seed{r(), r(), r(), r(), r(), r(), r(), r()};
        generator = mt19937_64(seed);
        unifx = uniform_real_distribution<>(bounds.getMinX(), bounds.getMaxX());
        unify = uniform_real_distribution<>(bounds.getMinY(), bounds.getMaxY());
    }

    Point generate() {
        return Point(unifx(generator), unify(generator));
    }
};

class alignas(128) Count {
public:
    long in = 0;
    long out = 0;
};

int main(int argc, char ** argv) {
    long N = 5E8;
    double R = 2;
    Circle c1(0, 0, 2);
    Circle c2(0, 1, 2);
    Circle c3(1, 0, 2);
    Circle c4(1, 1, 2);
    vector<Circle> circles{c1, c2, c3, c4};
    Bounds bounds = Bounds::find(circles);
    int numThreads = 2;
    N *= numThreads;
    cout << "starting bounds: " << bounds.getMinX() << " " << bounds.getMaxX()<< " " << bounds.getMinY() << " " << bounds.getMaxY() << "\n";
    for(int i = 0; i < 2; ++i) {
        long in = 0;
        long out = 0;
        vector<thread> threads;
        vector<Count> counts(numThreads);
        vector<Bounds> empericalBounds(numThreads);
        for (int t = 0; t < numThreads; ++t) {
            threads.emplace_back([&, t]() {
                Count &count = counts[t];
                Generator generator(bounds);
                Bounds & bounds = empericalBounds[t];
                for (long i = 0; i < N / numThreads; ++i) {
                    Point p = generator.generate();
                    bool insideAll = Circle::insideAll(p, circles);
                    if (insideAll) {
                        bounds.observed(p);
                        ++count.in;
                    } else {
                        ++count.out;
                    }
                }
                cout << "Observed bounds x: " << bounds.getMinX() << " " << bounds.getMaxX() << " y: " << bounds.getMinY() << " " << bounds.getMaxY() << "\n";
            });
        }
        for (thread &t : threads) {
            t.join();
        }
        for (int t = 0; t < numThreads; ++t) {
            in += counts[t].in;
            out += counts[t].out;
        }
        cout << N << " points: " << in / (double) (in + out) * bounds.getVolume() << " volume: " << R * R * M_PI << "\n";
        bounds = empericalBounds[0].grow(1.01);
    }
    return 0;
}