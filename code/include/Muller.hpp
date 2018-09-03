/**
 * Copyright (C) 2017-2018
 * Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 *
 * @author Pablo Alvarado
 * @date   18.08.2018
 */

#ifndef ANPI_MULLER_HPP
#define ANPI_MULLER_HPP

#include <vector>
#include <type_traits>
#include <math.h>

#include "PolynomialFormulaFormat.hpp"
#include <Deflation.hpp>

#include <boost/type_traits/is_complex.hpp>
#include <boost/math/tools/polynomial.hpp>

namespace anpi
{

/////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////

/// Enumerator makes explicit polishing roots
enum PolishEnum
{
  DoNotPolish,
  PolishRoots
};

namespace bmt = boost::math::tools; // for polynomial

//Muller method for one root
/**
   * Compute one root of the given polynomial using the Muller method.
   * All computations are done in complex numbers so we can get complex roots
   * result is a complex number with type of the inner type of U whis is double or float
   * @param[in] poly polynomial to be analyzed for roots
   * @param[out] result the root found
   * @param[in] start initial point for finding the root
   *
   */
template <class T, class U>
void muller(const bmt::polynomial<T> &poly,
            std::complex<typename anpi::detail::inner_type<U>::type> &result,
            const U start = U())
{
  //typedef to simplify instantiations of complex numers
  typedef typename anpi::detail::inner_type<U>::type type;
  typedef typename std::complex<type> complex;

  //maximun allowed iterations
  const unsigned int MAX_ITERATIONS = 50;
  const type eps = std::numeric_limits<type>::epsilon();
  //constant complex that are added
  const complex comp1 = complex(1), comp2 = complex(2), comp4 = complex(4);

  //simple Muller

  //choosing 3 equally distant points to start approximation
  //the choice to make the spacing 1 is arbitrary
  complex x1 = complex(start), x2 = x1 + comp1, x3 = x1 + comp2;

  //variables to hold calculations
  complex q, A, B, C, xi;
  //variables to hold the evaluated polynomials
  complex f1, f2, f3;

  size_t psize = poly.size();

  //size of the polynomial is 0
  if (psize < 2)
    return;

  for (unsigned int i = 0; i < MAX_ITERATIONS; i++)
  {

    //evaluate polynomials at the three given approximations
    f1 = bmt::evaluate_polynomial(&poly.data()[0], x1, psize);
    f2 = bmt::evaluate_polynomial(&poly.data()[0], x2, psize);
    f3 = bmt::evaluate_polynomial(&poly.data()[0], x3, psize);

    //after evaluating the functions we check if we have arrived at a 0
    if (std::abs(f1) <= eps)
    { //*****************revisar, creo q esta incorrecto
      // return x1;
      result = x1;
      return;
    }

    //calculating Mullers formulas
    q = (x1 - x2) / (x1 - x2);
    A = q * f1 - q * (comp1 + q) * f2 + q * q * f3;
    B = (comp2 * q + comp1) * f1 - (comp1 + q) * (comp1 + q) * f2 + q * q * f3;
    C = (comp1 + q) * f1;

    //calculate new approximation
    xi = B.real() >= 0 ? x1 - (x1 - x2) * (comp2 * C) / (B + std::sqrt(B * B - comp4 * A * C))
                       : x1 - (x1 - x2) * (comp2 * C) / (B - std::sqrt(B * B - comp4 * A * C));

    //set values for next iteration
    x3 = x2;
    x2 = x1;
    x1 = xi;

  } //end for
  //return

  throw("Mullers method for one Root::No root found");
} //end Muller()

/**
   * Compute the roots of the given polynomial using the Muller method.
   * @param[in] poly polynomial to be analyzed for roots
   * @param[out] roots all roots found
   * @param[in] start initial point for finding the roots
   * @param[in] polish indicate if polishing is needed or not.
   *
   * @return the number of roots found
   */
template <class T, class U>
void muller(const bmt::polynomial<T> &poly,
            std::vector<U> &roots,
            const PolishEnum polish = DoNotPolish,
            const U start = U())
{

  static_assert(std::is_floating_point<T>::value ||
                    boost::is_complex<T>::value,
                "T must be floating point or complex");
  static_assert(std::is_floating_point<U>::value ||
                    boost::is_complex<U>::value,
                "U must be floating point or complex");

  //we have to cast the inner type of midresult to the Inner type of T
  typedef typename anpi::detail::inner_type<T>::type Ttype;

  //define the polinomial to deflate to get the subsequent roots
  //and the polinomial to use for division
  typename bmt::polynomial<T> dfltPoly = poly, resDef2;
  T residuoPoly;
  //variable to store the result of the Muller method
  typedef typename anpi::detail::inner_type<U>::type type;

  std::complex<type> midResult;

  // the polynomial has to be of order 1 or above in order to be deflated
  //and for Muller to be able to find a root
  while (dfltPoly.size() > 1)
  {
    //calculate a root
    anpi::muller<T, U>(dfltPoly, midResult, start);

    //if polishing is on, polish that root with muller again
    if (polish == PolishRoots)
      anpi::muller<T, std::complex<type>>(poly, midResult, midResult);

    // si hay un NaN en los resultados
    if (midResult.real() == NAN)
      return;
    // ***************************

    //if Muller's method result was real
    if (anpi::detail::is_real(midResult))
    {
      //if the coefficients of the polynomial are also real
      if (std::is_floating_point<T>::value)
      {
        //we have to cast the result to the same values as the coefficients for deflate to work
        anpi::deflate<T>(dfltPoly, T(midResult.real()), residuoPoly);
      }

      //the coefficients are complex
      if (boost::is_complex<T>::value) // T is complex
      {
        anpi::deflate<T>(dfltPoly, Ttype(midResult.real()), residuoPoly);
      }

      //we are looking for real roots
      if (std::is_floating_point<U>::value)
      {
        //add to the results
        roots.push_back(midResult.real());
      }
      //looking for complex and real roots
      if (boost::is_complex<U>::value)
      {
        roots.push_back(midResult.real());
      }
    }
    else //Muller Result is complex
    {

      //if the coefficients of the polynomial are real
      if (std::is_floating_point<T>::value)
      {
        //special case of deflation, conjugate complex roots

        anpi::deflate2<T>(dfltPoly, std::complex<Ttype>(Ttype(midResult.real()), Ttype(midResult.imag())), resDef2);
      }
      if (boost::is_complex<T>::value)
      {
        //the coefficients are complex
        anpi::deflate<T>(dfltPoly, Ttype(midResult.real()), residuoPoly);
      }

      //if we are looking for complex and real roots
      if (boost::is_complex<U>::value)
      {
        roots.push_back(midResult.real());
      }
    } //end root is complex
  }   //end while
} //end Muller()

} // namespace anpi

#endif
