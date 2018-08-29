/**
 * Copyright (C) 2017-2018
 * Área Académica de Ingeniería en Computadoras, TEC, Costa Rica
 *
 * This file is part of the CE3102 Numerical Analysis lecture at TEC
 */


#include <boost/test/unit_test.hpp>

#include <iostream>
#include <exception>
#include <cstdlib>
#include <complex>
#include <string>
#include <sstream>

#include <PolynomialFormulaFormat.hpp>

namespace bmt=boost::math::tools; // for polynomial

// normal allocator

typedef std::complex<double> dcomplex;
typedef std::complex<float>  fcomplex;

BOOST_AUTO_TEST_SUITE( PolynomialFormulaFormat )

BOOST_AUTO_TEST_CASE( PolyFormulaFormat ) {
  bmt::polynomial<float> p1 = {{1.f,0.f,-2.f,0.f,0.f,3.f}};
  std::string str = "3x^5 - 2x^2 + 1";
  BOOST_CHECK(anpi::polynomialFormulaFormat(p1)==str);

  typedef std::complex<double> dcomplex;

  bmt::polynomial< dcomplex > p2 =
    {{1.,0.,dcomplex(-2,1),0.,0.,dcomplex(0,3.f)}};
  str = "(0,3)x^5 + (-2,1)x^2 + 1";
  BOOST_CHECK(anpi::polynomialFormulaFormat(p2)==str);
}
  
BOOST_AUTO_TEST_SUITE_END()
