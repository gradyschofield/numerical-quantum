//
// Created by Grady Schofield on 4/23/21.
//

#include<random>
#include<thread>
#include<iostream>

using namespace std;

class Point {
    double x, y, z;
public:
    Point(double x, double y, double z) : x(x), y(y), z(z){}

    double getX() const {
        return x;
    }

    double getY() const {
        return y;
    }

    double getZ() const {
        return z;
    }
};

class Sphere {
    double x, y, z, r, r2;
public:
    Sphere(double x, double y, double z, double r)
            : x(x), y(y), z(z), r(r), r2(r*r)
    {
    }

    bool inside(Point const & p) {
        double dx = p.getX() - Sphere::x;
        double dy = p.getY() - Sphere::y;
        double dz = p.getZ() - Sphere::z;
        return dx*dx + dy*dy + dz*dz < r2;
    }

    static bool insideAll(Point const & p, vector<Sphere> const & spheres) {
        for(Sphere const & s : spheres) {
            double dx = p.getX() - s.getX();
            double dy = p.getY() - s.getY();
            double dz = p.getZ() - s.getZ();
            if(dx*dx + dy*dy + dz*dz> s.getR2()) {
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

    double getZ() const {
        return z;
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
    double minz = numeric_limits<double>::max();
    double maxz = -numeric_limits<double>::max();

public:
    Bounds(){
    }

    Bounds(double minx, double maxx,
           double miny, double maxy,
           double minz, double maxz)
            : minx(minx), maxx(maxx), miny(miny), maxy(maxy), minz(minz), maxz(maxz)
    {
    }

    void observed(Point const & p) {
        minx = min(minx, p.getX());
        maxx = max(maxx, p.getX());
        miny = min(miny, p.getY());
        maxy = max(maxy, p.getY());
        minz = min(minz, p.getZ());
        maxz = max(maxz, p.getZ());
    }

    static Bounds find(vector<Sphere> const & spheres) {
        double minx, miny, minz = miny = minx = numeric_limits<double>::max();
        double maxx, maxy, maxz = maxy = maxx = -numeric_limits<double>::max();
        for(Sphere const & c : spheres) {
            minx = min(minx, c.getX() - c.getR());
            maxx = max(maxx, c.getX() + c.getR());
            miny = min(miny, c.getY() - c.getR());
            maxy = max(maxy, c.getY() + c.getR());
            minz = min(minz, c.getZ() - c.getR());
            maxz = max(maxz, c.getZ() + c.getR());
        }
        return Bounds(minx, maxx, miny, maxy, minz, maxz);
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

    double getMinZ() const {
        return minz;
    }

    double getMaxZ() const {
        return maxz;
    }

    double getVolume() const {
        return (maxx - minx) * (maxy - miny) * (maxz - minz);
    }

    Bounds grow(double p) {
        double midx = 0.5*(maxx + minx);
        double midy = 0.5*(maxy + miny);
        double midz = 0.5*(maxz + minz);
        double widthx = 0.5*(maxx - minx);
        double widthy = 0.5*(maxy - miny);
        double widthz = 0.5*(maxz - minz);
        return Bounds(midx - widthx * p,
                      midx + widthx * p,
                      midy - widthy * p,
                      midy + widthy * p,
                      midz - widthz * p,
                      midz + widthz * p);
    }
};

class Generator {
    mt19937_64 generator;
    uniform_real_distribution<> unifx;
    uniform_real_distribution<> unify;
    uniform_real_distribution<> unifz;
public:
    Generator(Bounds const & bounds) {
        random_device r;
        seed_seq seed{r(), r(), r(), r(), r(), r(), r(), r()};
        generator = mt19937_64(seed);
        unifx = uniform_real_distribution<>(bounds.getMinX(), bounds.getMaxX());
        unify = uniform_real_distribution<>(bounds.getMinY(), bounds.getMaxY());
        unifz = uniform_real_distribution<>(bounds.getMinZ(), bounds.getMaxZ());
    }

    Point generate() {
        return Point(unifx(generator), unify(generator), unifz(generator));
    }
};

class alignas(128) Count {
public:
    long in = 0;
    long out = 0;
};

int main(int argc, char ** argv) {
    long N = 1E6;
    double R = 2;
    Sphere s1(0, 0, 0, 2);
    Sphere s2(0, 1, 0.5, 2);
    Sphere s3(1, 0, 0, 2);
    Sphere s4(1, 1, 0, 2);
    vector<Sphere> spheres{s1, s2, s3, s4};
    Bounds bounds = Bounds::find(spheres);
    int numThreads = 2;
    N *= numThreads;
    cout << "starting bounds: " << bounds.getMinX() << " " << bounds.getMaxX()
        << " " << bounds.getMinY() << " " << bounds.getMaxY()
        << " " << bounds.getMinZ() << " " << bounds.getMaxZ() << "\n";
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
                    bool insideAll = Sphere::insideAll(p, spheres);
                    if (insideAll) {
                        bounds.observed(p);
                        ++count.in;
                    } else {
                        ++count.out;
                    }
                }
                cout << "Observed bounds x: " << bounds.getMinX() << " " << bounds.getMaxX()
                    << " y: " << bounds.getMinY() << " " << bounds.getMaxY()
                    << " z: " << bounds.getMinZ() << " " << bounds.getMaxZ() << "\n";
            });
        }
        for (thread &t : threads) {
            t.join();
        }
        for (int t = 0; t < numThreads; ++t) {
            in += counts[t].in;
            out += counts[t].out;
        }
        cout << N << " points: " << in / (double) (in + out) * bounds.getVolume() << " volume: " << 4 * R * R * R * M_PI / 3 << "\n";
        bounds = empericalBounds[0].grow(1.05);
        N = 5E7;
    }
    return 0;
}
