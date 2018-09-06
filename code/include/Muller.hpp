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

#include <iostream>
#include <vector>
#include <type_traits>
#include <math.h>
#include <complex.h>

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

////////////////////////////////////////////////////////////////////////////////////////////

template <class T>
inline void
addResult(std::vector<std::complex<T>> &results, std::complex<T> &res)
{
  results.push_back(res);
}

template <class T>
inline void
addResult(std::vector<T> &results, std::complex<T> &res)
{
  results.push_back(res.real());
}
/////////////////////////////////////////////////////////////////////////////////////////////////
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
  typedef typename anpi::detail::inner_type<U>::type Utype;
  typedef typename std::complex<Utype> complex;

  //maximun allowed iterations
  const unsigned int MAX_ITERATIONS = 100000;
  const Utype eps = std::numeric_limits<Utype>::epsilon();
  //estimated error
  // complex ea;
  //constant complex numbers that are added in the equations
  const complex comp1 = complex(1), comp2 = complex(2), comp4 = complex(4), compZero(0);

  //choosing 3 equally distant points to start approximation
  //the choice to make the spacing 1 is arbitrary
  complex x1 = complex(start), x2 = x1 + comp1, x3 = x1 + comp2;

  //variables to hold calculations
  complex q, A, B, C, plusDenom, minusDenom, plusxi, minusxi, xi;
  //variables to hold the evaluated polynomials
  complex f1, f2, f3;

  //size and degree of the polynomial
  size_t psize = poly.size();
  int deg = anpi::polyDegree(poly);

  //degree of the polynomial is 0, there are no root
  //for this is a constant
  if (deg < 1)
  {
    result = std::numeric_limits<T>::quiet_NaN();
    return;
  }

  //first evaluation
  //evaluate polynomials at the given approximations
  f2 = bmt::evaluate_polynomial(&poly.data()[0], x2, psize);
  f3 = bmt::evaluate_polynomial(&poly.data()[0], x3, psize);

  for (unsigned int i = 0; i < MAX_ITERATIONS; ++i)
  {

    //evaluate polynomials at the three given approximations
    f1 = bmt::evaluate_polynomial(&poly.data()[0], x1, psize);

    //after evaluating the functions we check if we have arrived at a 0a
    if (std::abs(f1) <= eps)
    {
      //check if imaginary part is too small
      if (0 < x1.imag() && x1.imag() < 2 * eps)
      {
        x1.imag(0);
      }
      if (0 < x1.real() && x1.real() < 2 * eps)
      {
        x1.real(0);
      }
      //check if real part is too small

      // return x1;

      result = x1;
      std::cout << "Found on " << i << " iterations" << std::endl;
      return;
    }

    //calculating Mullers formulas
    q = (x1 - x2) / (x2 - x3);
    A = q * f1 - q * (comp1 + q) * f2 + q * q * f3;
    B = (comp2 * q + comp1) * f1 - (comp1 + q) * (comp1 + q) * f2 + q * q * f3;
    C = (comp1 + q) * f1;

    plusDenom = B + std::sqrt((B * B) - (comp4 * A * C));
    plusxi = x1 - ((x1 - x2) * (comp2 * C) / plusDenom);
    minusDenom = B - std::sqrt((B * B) + (comp4 * A * C));
    minusxi = x1 - ((x1 - x2) * (comp2 * C) / minusDenom);
    //calculate new approximation
    xi = (std::abs(plusDenom) > std::abs(minusDenom) && std::abs(x1 - plusxi) < std::abs(x1 - minusxi))
             ? plusxi
             : minusxi;

    //our estimation became a NaN value, no need to continue computing
    if (std::isnan(std::abs(xi)))
    {
      std::cout << "xi se hizo NaN " << std::endl;
      result = std::numeric_limits<T>::quiet_NaN();
      return;
    }
    //set values for next iteration
    x3 = x2;
    x2 = x1;
    x1 = xi;

    f3 = f2;
    f2 = f1;

  } //end for
  //return

  std::cout << "no se encontro en Max Iterations = " << MAX_ITERATIONS << std::endl;
  //no Root found on MAX_ITERATIONS
  result = std::numeric_limits<T>::quiet_NaN();

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

  const bool isRealCoeff = std::is_floating_point<T>::value;
  const bool isComplexRoots = boost::is_complex<U>::value;

  typedef typename anpi::detail::inner_type<T>::type Ttype;
  typedef typename anpi::detail::inner_type<U>::type Utype;

  //define the polinomial to deflate to get the subsequent roots
  //and the polinomial to use for division
  typename bmt::polynomial<T> dfltPoly = poly, resDef2;
  T residuoPoly;

  int countRoots = 0;

  //variable to hold the root every iteration
  std::complex<Utype> midResult;

  //degree of the polynomial
  int deg = anpi::polyDegree(poly);

  //polynomial order is less than 1 (a constant), no possible results
  if (deg < 1)
  {
    roots.push_back(std::numeric_limits<U>::quiet_NaN());
    return;
  }

  // the polynomial has to be of order 1 or above in order to be deflated
  //and for Muller to be able to find a root
  while (deg > 0)
  {
    //calculate a root
    anpi::muller<T, U>(dfltPoly, midResult, start);

    //if polishing is on, polish that root with muller again
    if (polish == PolishRoots)
      anpi::muller<T, std::complex<Utype>>(poly, midResult, midResult);

    if (isnan(abs(midResult)))
    {
      ++countRoots;
      roots.push_back(std::numeric_limits<U>::quiet_NaN());
      return;
    }

    //if Muller's method result was real
    if (anpi::detail::is_real(midResult))
    {
      //if the coefficients of the polynomial are also real
      if (isRealCoeff)
      {
        //we have to cast the result to the same values as the coefficients for deflate to work
        dfltPoly = anpi::deflate<T>(dfltPoly, T(midResult.real()), residuoPoly);
      }
      else // T is complex
      {

        dfltPoly = anpi::deflate<T>(dfltPoly, T(midResult.real()), residuoPoly);
      }

      //we are looking for real roots
      if (!isComplexRoots)
      {
        ++countRoots;
        //add to the results
        anpi::addResult(roots, midResult);
      }
      //looking for complex and real roots
      else //U is complex, roots<U> should be complex
      {
        ++countRoots;
        anpi::addResult(roots, midResult);
      }
    }
    else //Muller Result is complex
    {
      //if the coefficients of the polynomial are real
      if (isRealCoeff)
      {
        // //special case of deflation, conjugate complex roots
        // static_assert(std::is_floating_point<T>::value,
        //               "T must be floating point or complex");//static assert is failing

        dfltPoly = anpi::deflate2<T>(dfltPoly, std::complex<Ttype>(Ttype(midResult.real()), Ttype(midResult.imag())), resDef2);

        if (isComplexRoots)
        {
          ++countRoots;
          std::complex<Utype> conjugate = std::conj(midResult);
          anpi::addResult(roots, conjugate);
        }

        //if we do special deflation, we are removing 2 roots so substract an extra degree
        --deg;
      }
      else
      {
        //the coefficients are complex
        dfltPoly = anpi::deflate<T>(dfltPoly, Ttype(midResult.real()), residuoPoly);
      }

      //if we are looking for complex and real roots
      if (isComplexRoots)
      {
        ++countRoots;
        // complexResults.push_back(midResult);
        anpi::addResult(roots, midResult);
      }
    } //end root is complex

    //in order not to calculate the degree again we just substract one because we deflated the polynomial
    --deg;
  } //end while

  if (countRoots == 0)
  {
    roots.push_back(std::numeric_limits<U>::quiet_NaN());
    return;
  }
  // if (std::is_same<U, std::complex<double>>::value || std::is_same<U, std::complex<float>>::value)
  // {
  //   roots = complexResults;
  // }
  // else
  // roots = realResults;

} //end Muller()

/////////////////////////////////////////////////////////////

} // namespace anpi

#endif
