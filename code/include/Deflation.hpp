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


namespace anpi {
  namespace bmt=boost::math::tools; // bt as alias for boost::math::tools
  
  /**
   * Deflate polynomial
   *
   * @param[in] poly Input polynomial to be deflated with the provided root
   * @param[in] root Root of poly to be deflated with.
   * @param[out] residuo Residual of the polynomial deflation
   * @return deflated polynomial
   */
  template<class T>
  bmt::polynomial<T> deflate(const bmt::polynomial<T>& poly,
                             const T& root,
                             T& residuo,
                             T& tolerance=anpi::epsilon<T>());

  /**
   * Deflate polynomial with a second order polynomial.
   *
   * The second order polynomial equals x^2 -2 Re(root)x + |root|^2.
   *
   * @param[in] poly Input polynomial to be deflated with the provided root
   * @param[in] root Root of poly to be deflated with.
   * @param[out] residuo Residual of the polynomial deflation
   * @return deflated polynomial
   */
  template<class T>
  bmt::polynomial<T> deflate2(const bt::polynomial<T>& poly,
                              const std::complex<T>& root,
                              bt::polynomial<T>& residuo,
                              T& tolerance=anpi::epsilon<T>());

  
}



#endif
