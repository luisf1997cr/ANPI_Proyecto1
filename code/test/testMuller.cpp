#include <vector>
#include <iostream>
#include <iomanip>
#include <complex>

// #include "PolynomialFormulaFormat.hpp"
#include <boost/type_traits/is_complex.hpp>
#include <boost/math/tools/polynomial.hpp>
#include <cmath>
#include <type_traits>

#include <Muller.hpp>
#include <JenkinsTraub.hpp>

namespace bmt = boost::math::tools;

int main()
{
    //Variables
    std::complex<float> cf0 = 0.f, cf1 = 1.f, cf2 = 2.f, cf3 = 3.f, cf4 = 4.f, rescf, resMCF;
    std::complex<double> cd0 = 0., cd1 = 1., cd2 = 2., cd3 = 3., cd4 = 4., rescd, resMCD;
    double d0 = 0, d1 = 1, d2 = 2, d3 = 3, d4 = 4, resd, resMD;
    float f0 = 0, f1 = 1, f2 = 2, f3 = 3, f4 = 4, resf, resMF;

    //___Polinomio poly1 = 2X^2+4X_________roots: 0, -2_________________________________________________________
    bmt::polynomial<std::complex<float>> poly1cf = {cf0, cf4, cf2};
    bmt::polynomial<std::complex<double>> poly1cd = {cd0, cd4, cd2};
    bmt::polynomial<double> poly1d = {d0, d4, d2};
    bmt::polynomial<float> poly1f = {f0, f2, f4};

    anpi::deflate<std::complex<float>>(poly1cf, -cf2, rescf);
    anpi::deflate<std::complex<double>>(poly1cd, -cd2, rescd);
    anpi::deflate<double>(poly1d, -d2, resd);
    anpi::deflate<float>(poly1f, -f2, resf);

    anpi::muller<std::complex<float>, std::complex<float>>(poly1cf, resMCF, cf1);
    anpi::muller<std::complex<double>, std::complex<double>>(poly1cd, resMCD, cd1);
    anpi::muller<double, double>(poly1d, resMD, d1);
    anpi::muller<float, float>(poly1f, resMF, f1);

    //____x^3-x^2+4x-4___roots = 1, +-2i_________________________________________________________________________
    bmt::polynomial<std::complex<float>> poly2cf = {-cf4, cf4, -cf1, cf1};
    bmt::polynomial<std::complex<double>> poly2cd = {-cd4, cd4, -cd1, cd1};
    bmt::polynomial<double> poly2d = {-d4, d4, -d1, d1};
    bmt::polynomial<float> poly2f = {-f4, f4, -f1, f1};

    anpi::deflate<std::complex<float>>(poly2cf, cf1, rescf);
    anpi::deflate<std::complex<double>>(poly2cd, cd1, rescd);
    anpi::deflate<double>(poly2d, d1, resd);
    anpi::deflate<float>(poly2f, f1, resf);

    anpi::muller<std::complex<float>, std::complex<float>>(poly2cf, resMCF, cf1);
    anpi::muller<std::complex<double>, std::complex<double>>(poly2cd, resMCD, cd1);
    // anpi::muller<double, double>(poly2d, resMD, d1);
    // anpi::muller<float, float>(poly2f, resMF, f1);

    //______x⁴-7x³+13x²+23x-78___roots = -2,3,3+2i,3-2i_____
    bmt::polynomial<std::complex<float>> poly3cf = {78.f, 23.f, 13.f, -7.f, cf1};
    // bmt::polynomial<std::complex<double>> poly3cd = {78d, 23d, 13d, -7d, cd1};
    bmt::polynomial<double> poly3d = {78, 23, 13, -7, d1};
    bmt::polynomial<float> poly3f = {78.f, 23.f, 13.f, -7.f, f1};

    anpi::deflate<std::complex<float>>(poly3cf, -cf2, rescf);
    // anpi::deflate<std::complex<double>>(poly3cd, -cd2, rescf);
    anpi::deflate<double>(poly3d, -d2, resd);
    anpi::deflate<float>(poly3f, -f2, resf);

    anpi::muller<std::complex<float>, std::complex<float>>(poly3cf, resMCF, cf1);
    // anpi::muller<std::complex<double>, std::complex<double>>(poly3cd, resMCD, cd1);
    // anpi::muller<double, double>(poly3d, resMD, d1);
    // anpi::muller<float, float>(poly3f, resMF, f1);

    return 0;
}