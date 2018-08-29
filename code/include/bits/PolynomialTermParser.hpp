/**
 * Copyright (C) 2017-2018
 * Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 *
 * @author Pablo Alvarado
 * @date   17.08.2018
 */

#ifndef ANPI_POLYNOMIAL_TERM_PARSER_HPP
#define ANPI_POLYNOMIAL_TERM_PARSER_HPP

#include <iostream>
#include <string>
#include <complex>

#include <boost/config/warning_disable.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_object.hpp>
#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/fusion/include/io.hpp>
#include <boost/bind.hpp>
#include <boost/spirit/include/phoenix.hpp>

#include <boost/math/tools/polynomial.hpp>

namespace anpi {
  namespace detail {
    
    namespace qi = boost::spirit::qi;
    namespace ascii = boost::spirit::ascii;

    /**
     * Polynomial term coefficient x^exponent
     *
     * The most general case handles complex coefficients and integral
     * exponents.
     */
    struct Term {

      Term() : coefficient(1.0), xExponent(0) {}
      
      /**
       * Coefficient of the polynomial term.
       *
       * If not given, it is assumed as 1
       */
      std::complex<double> coefficient;

      /** 
       * Exponent of the polynomial term.
       *
       * If not given, it is assumed as 0
       */
      unsigned int xExponent;

      /**
       * Compare if this term is less than the other term based on their
       * exponents only
       */
      inline bool operator<(const Term& other) const {
	return xExponent < other.xExponent;
      }

      /**
       * Compare if this term is less than the other term based on their
       * exponents only
       */
      inline bool operator<=(const Term& other) const {
	return xExponent <= other.xExponent;
      }

    };

    /**
     * The parser for polynomial expresions
     */
    template <typename Iterator>
    struct PolynomialTermParser : qi::grammar<
                                           Iterator,
                                           std::vector<Term>(),
                                           ascii::space_type
                                          >
    {
      /// Construct the parser 
      PolynomialTermParser() : PolynomialTermParser::base_type(start) {
        namespace phx = boost::phoenix;
	
        using qi::uint_;
        using qi::double_;
        using qi::_val;
        using qi::eps;
	
	// The exponential part is restricted to integrals after ^
        exponent %= '^' >> uint_;

	// An x term can optionally have an exponent appended
        x_term %= 'x' >> ( exponent | eps[_val = 1] );

	// sqrt(-1)
	static const std::complex<double> i(0,1);
	
	// Parse a complex number as (a,b) (a real part, b imaginary part)
	// accept also a and (a)
        cplx %= '(' >> double_
		    >> -(','
		    >> double_[_val+=qi::_1*i])
		    >> ')'
	  | double_;
        

	// polynomial term has an optional coefficient and an optional term
	poly_term =
	  eps[_val = Term()]
          >> -cplx[phx::bind(&Term::coefficient, _val) = qi::_1]
          >> -x_term[phx::bind(&Term::xExponent, _val) = qi::_1];
	
	// a negated term like before but, uhm, negated
        negative_poly_term =
	  eps[_val = Term()]
	  >>
          (
           cplx[phx::bind(&Term::coefficient, _val) =  -1.0 * qi::_1]
           | eps[phx::bind(&Term::coefficient, _val) = -1.0]
          )
	  >>
          -x_term[phx::bind(&Term::xExponent, _val) = qi::_1];
	
	// the final polynomial as a sequence of terms
        start = eps[_val = std::vector<Term>()]
          >> poly_term[phx::push_back(_val, qi::_1)]
          >> *(
               ('+' >> poly_term[phx::push_back(_val, qi::_1)])
               |
               ('-' >> negative_poly_term[phx::push_back(_val, qi::_1)])
              );
      }

      // All rules
      qi::rule<Iterator, unsigned int(), ascii::space_type> exponent;
      qi::rule<Iterator, unsigned int(), ascii::space_type> x_term;
      qi::rule<Iterator, std::complex<double>(), ascii::space_type> cplx;
      qi::rule<Iterator, Term(), ascii::space_type> poly_term;
      qi::rule<Iterator, Term(), ascii::space_type> negative_poly_term;
      qi::rule<Iterator, std::vector<Term>(), ascii::space_type> start;
    };
  }
}

#endif
