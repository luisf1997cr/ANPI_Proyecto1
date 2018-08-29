/**
 * Copyright (C) 2017-2018
 * Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 *
 * @author Pablo Alvarado
 * @date   17.08.2018
 */

#ifndef ANPI_POLYNOM_PARSER_HPP
#define ANPI_POLYNOM_PARSER_HPP

#include <iostream>
#include <string>

#include <boost/math/tools/polynomial.hpp>

#include "bits/PolynomialTermParser.hpp"
#include "Exception.hpp"

namespace anpi {

  namespace bmt=boost::math::tools; // for polynomial

  namespace aux {
    // Auxiliary method to add the coefficients
    template<class T>
    inline typename std::enable_if<std::is_floating_point<T>::value,T>::type
    add(const T val,const std::complex<double>& coef) {
      return val+coef.real();
    }
    
    template<class T>
    inline std::complex<T> 
    add(const std::complex<T>& val,const std::complex<double>& coef) {
      return std::complex<T>(static_cast<T>(val.real()+coef.real()),
                             static_cast<T>(val.imag()+coef.imag()));
    }
  }
  
  /**
   * Parse a polynomial string and convert it to a boost polynomial.
   *
   * The variable of the polynomial string must be 'x'.
   *
   * Valid strings are:
   *
   * - x^2 + 4.5x - 3
   * - 3x^4 + x -2x^2 + 5
   * - (3,2)x^2 + (0,1)x^3 - 5x^4 + (2,0)
   *
   * where the (a,b) terms refer to complex numbers with real part a and 
   * imaginary part b.
   *
   * Errors will be reported via exceptions.
   */
  template<class T>
  bmt::polynomial<T> parsePolynomial(const std::string& str) {
    using boost::spirit::ascii::space;
    
    typedef std::string::const_iterator iterator_type;
    typedef anpi::detail::PolynomialTermParser<iterator_type> polynomial_parser;

    polynomial_parser pparser;

    // store the parsed terms here
    std::vector<anpi::detail::Term> poly;

    iterator_type start = str.begin();
    iterator_type end   = str.end();
    bool r = boost::spirit::qi::phrase_parse(start,end,pparser,space,poly);
    if (!r) {
      throw Exception("Could not parse given string");
    } else if (start != end) {
      throw Exception(std::string("Syntactic error at position ") +
		      boost::lexical_cast<std::string>(start-str.begin()));
    }

    // sort all terms by their exponent, so we can catch if there are
    // several terms with the same order or if there are missing exponents
    std::sort(poly.begin(),poly.end());
    unsigned int order = 1u + poly.back().xExponent;
    std::vector<T> allCoefs(order,T());
    
    anpi::detail::Term term;

    // The polynomial class expects the constant term first, and then
    // the coefficients of x,x^2,x^3 and so on
    for (auto& it : poly) {
      // extract the coefficient 
      const unsigned int xExp = it.xExponent;
      const std::complex<double> coef = it.coefficient;
      //
      if (std::is_floating_point<T>::value &&
	  (std::abs(coef.imag())>std::numeric_limits<double>::epsilon())) {
	throw Exception("Complex coefficients not allowed "
			"for real polynomials");
      } else {
	allCoefs[xExp]=aux::add(allCoefs[xExp],coef);
      }
    }

    return bmt::polynomial<T>(allCoefs.begin(),allCoefs.end());
  }
  
}


#endif
