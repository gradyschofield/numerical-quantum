//
// Created by grady on 6/22/20.
//

#ifndef SAMPLES_COMPLEX_H
#define SAMPLES_COMPLEX_H

#include<ostream>
#include<cmath>

template<typename T>
class Complex {
    T re;
    T im;

public:
    Complex()
        : re(0), im(0)
    {
    }

    Complex(T re)
        : re(re), im(0)
    {
    }

    Complex(T re, T im)
        : re(re), im(im)
    {
    }

    T getRe() const {
        return re;
    }

    T getIm() const {
        return im;
    }
};

template<typename T>
T abs(Complex<T> const & c1) {
    return sqrt(c1.getRe() * c1.getRe() + c1.getIm() * c1.getIm());
}

template<typename T>
T abs2(Complex<T> const & c1) {
    return c1.getRe() * c1.getRe() + c1.getIm() * c1.getIm();
}

template<typename T>
Complex<T> conjugate(Complex<T> const & c) {
    return Complex<T>(c.getRe(), -c.getIm());
}

template<typename T>
Complex<T> operator+(Complex<T> const & c1, Complex<T> const & c2) {
    return Complex<T>(c1.getRe() + c2.getRe(), c1.getIm() + c2.getIm());
}

template<typename T>
Complex<T> operator-(Complex<T> const & c1, Complex<T> const & c2) {
    return Complex<T>(c1.getRe() - c2.getRe(), c1.getIm() - c2.getIm());
}

template<typename T>
Complex<T> operator*(Complex<T> const & c1, Complex<T> const & c2) {
    return Complex<T>(c1.getRe() * c2.getRe() - c1.getIm() * c2.getIm(), c1.getRe() * c2.getIm() + c1.getIm() * c2.getRe());
}

template<typename T>
Complex<T> operator/(Complex<T> const & c1, Complex<T> const & c2) {
    return (c1 * conjugate(c2)) / abs2(c2);
}

template<typename T>
Complex<T> operator/(Complex<T> const & c, T r) {
    T rinv = 1.0 / r;
    return Complex<T>(c.getRe() * rinv, c.getIm() * rinv);
}

template<typename T>
Complex<T> exp(Complex<T> const & x) {
    T a = exp(x.getRe());
    T c = cos(x.getIm());
    T s = sin(x.getIm());
    return Complex<T>(a*c, a*s);
}

template<typename T>
std::ostream & operator<<(std::ostream & os, Complex<T> const & c) {
    os << c.getRe() << (c.getIm() < 0 ? "-" : "+");
    if(fabs(c.getIm()) != 1) {
        os << fabs(c.getIm());
    }
    os << "i";
    return os;
}

#endif //SAMPLES_COMPLEX_H
