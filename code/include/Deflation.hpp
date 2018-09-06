/**
 * Copyright (C) 2017-2018
 * Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 *
 * @author Pablo Alvarado
 * @date   18.08.2018
 */

#ifndef ANPI_DEFLATION_HPP
#define ANPI_DEFLATION_HPP

#include <vector>
#include <type_traits>

#include <boost/type_traits/is_complex.hpp>
#include <boost/math/tools/polynomial.hpp>

namespace anpi
{
namespace bmt = boost::math::tools; // bt as alias for boost::math::tools

/**Function to know what is the degree of the polynomial
  *meaning, what is the highest power of x
  *
  *
  *
  */

template <class T>
int polyDegree(const bmt::polynomial<T> &poly)
{
  int deg = 0;
  int psize = poly.size();
  for (int i = 0; i < psize; ++i)
  {
    if (poly[i] != T(0))
      deg = i;
  }

  return deg;
}

/**
   * Divide a polynomial u by a polynomial v, and return the quotient and remainder polynomials
   * q and r, respectively. The four polynomials are represented as vectors of coefficients, each
   * starting with the constant term. There is no restriction on the relative lengths of u and v, and
   * either may have trailing zeros (represent a lower degree polynomial than its length allows). q
   * and r are returned with the size of u, but will usually have trailing zeros.
   *
   *
   * @param[in] u polynomial to be divided
   * @param[in] v polynomial to divide by
   * @param[out] q quotient of the polynomial division
   * @param[out] r residual of the polynomial division
   * 
   */
template <class T>
void poldiv(bmt::polynomial<T> &u, bmt::polynomial<T> &v, bmt::polynomial<T> &q, bmt::polynomial<T> &r)

{
  int degU = anpi::polyDegree<T>(u), degV = anpi::polyDegree<T>(v);

  int k, j;
  while (degV >= 0 && v[degV] == T(0))
    --degV;
  if (degV < 0)
    throw("poldiv divide by zero polynomial");
  //initialize to set size same as the polynomial we are dividing
  q = u;
  r = u;

  q[degU] = T(0);

  for (k = degU - degV; k >= 0; --k)
  {
    q[k] = r[degV + k] / v[degV];
    for (j = degV + k - 1; j >= k; --j)
      r[j] -= q[k] * v[j - k];
  }
  for (j = degV; j <= degU; ++j)
    r[j] = T(0);
}

/**
   * Deflate polynomial
   *
   * @param[in] poly Input polynomial to be deflated with the provided root
   * @param[in] root Root of poly to be deflated with.
   * @param[out] residuo Residual of the polynomial deflation
   * @return deflated polynomial
   */
template <class T>
bmt::polynomial<T> deflate(const bmt::polynomial<T> &poly,
                           const T &root,
                           T &residuo)
{
  typename bmt::polynomial<T> result = poly;
  unsigned int deg = anpi::polyDegree(poly);
  residuo = poly[deg];
  result[deg] = T(0);
  T swap;
  for (int i = deg - 1; i >= 0; --i)
  {
    swap = result[i];
    result[i] = residuo;
    residuo = swap + residuo * root;
  }

  return result;
}

/**
   * Deflate polynomial with a second order polynomial. Coefficients of the
   * polynomial have to be Real
   *
   * The second order polynomial equals x^2 -2 Re(root)x + |root|^2.
   *
   * @param[in] poly polynomial to be deflated with the provided root
   * @param[in] root complex Root of polynomial to be deflated with.
   * @param[out] residuo Residual of the polynomial deflation
   * @return deflated polynomial
   */
template <class T>
bmt::polynomial<T> deflate2(const bmt::polynomial<T> &poly,
                            const std::complex<T> &root,
                            bmt::polynomial<T> &residuo)

{
  // create a polynomials with 1 as coefficients, of size the same as the
  // polynomial we want to divide to hold the quotient and reminder. Also
  //a denominator
  // T psize = static_cast<T>(poly.size());
  typename bmt::polynomial<T> qnt = poly, rmdr = poly, pol = poly, dnmtr;

  T a = root.real(),
    b = root.imag();
  dnmtr = {(a * a + b * b), -(T(2) * a), T(1)};

  anpi::poldiv(pol, dnmtr, qnt, rmdr);

  residuo = rmdr;

  return qnt;
}

} // namespace anpi

#endif
