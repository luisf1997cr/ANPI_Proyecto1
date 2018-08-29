/**
 * Copyright (C) 2017-2018
 * Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 *
 * @author Pablo Alvarado
 * @date   17.08.2018
 */

#ifndef ANPI_POLYNOMIAL_FORMULA_FORMAT_HPP
#define ANPI_POLYNOMIAL_FORMULA_FORMAT_HPP

#include <iostream>
#include <string>

#include <boost/math/tools/polynomial.hpp>
#include <boost/lexical_cast.hpp>

namespace anpi {

  namespace bmt=boost::math::tools; // for polynomial

  namespace detail {
    
    /**
     * Metafunction to extract the inner type defined as the value_type of
     * a complex type or the type itself otherwise.
     */
    template<class T>
    struct inner_type { // Default implementation returns the type itself
      typedef T type;
    };
    
    template<class T>
    struct inner_type< std::complex<T> > { // specialization for complex types
      typedef T type;
    };
    
    /**
     * Function to check if a number is real, even though it may be
     * a std::complex (i.e. with imag() == 0)
     */
    inline constexpr bool is_real(const double& val) { return true; }
    inline constexpr bool is_real(const float& val)  { return true; }
    
    template <typename T>
    inline bool is_real(const std::complex<T>& val)  {
      return val.imag() == T(0);
    }
    
    /**
     * Function to check if a number is negative
     *
     * A number is negative if it is real and negative.
     */
    template<typename T>
    bool is_negative(const T& val) {
      typedef typename inner_type<T>::type type;
      return is_real(val) && (std::real(val) < type(0));
    }
    
    
    
    template <typename T>
    std::string signStr(const T &x) {
      return is_negative(x) ? "-" : "+";
    }
    
    
    template<class T>
    std::string justCoefficient(const T& x) {
      std::string result;
      if (!is_real(x)) {
        result += boost::lexical_cast<std::string>(x);
      } else {
        if (is_negative(x)) {
          result += "-";
        }
        if (std::abs(x) != T(1)) {
          result += boost::lexical_cast<std::string>(std::abs(std::real(x)));
        }
      }
      
      return result;
    }

    enum ShowOnesEnum {
      DoNotShowOnes,
      ShowOnes
    };
    
    template<class T>
    std::string innerCoefficient(const T& x,
                                 ShowOnesEnum showOnes=DoNotShowOnes) {
      
      std::string result(" " + signStr(x) + " ");
      
      if (!is_real(x)) {
        result += boost::lexical_cast<std::string>(x);
      } else {
        if ( (showOnes==ShowOnes) || (std::abs(x) != T(1)) ) {
          result += boost::lexical_cast<std::string>(std::abs(std::real(x)));
        }
      }
      return result;
    }
  }
  
  /**
   * Output in formula format.
   *
   * For example: 
   * from a polynomial in Boost container storage  [ 10, -6, -4, 3 ]
   * show as human-friendly formula notation: 3x^3 - 4x^2 - 6x + 10.
   */
  template <typename T>
  std::string polynomialFormulaFormat(bmt::polynomial<T> const &a) {
    
    // Empty polynomial just shown as zero
    if (a.size() == 0) {
      return boost::lexical_cast<std::string>(T(0));
    }
    
    // First one is a special case as it may need unary negate.
    unsigned i = a.size() - 1;

    while ( (i>0) && (a[i] == T()) ) {
      --i;
    }
    std::string result = detail::justCoefficient(a[i]);

    
    if (i > 0)  {
      result += "x";
      if (i > 1) {
        result += "^" + boost::lexical_cast<std::string>(i);
        i--;
        for (; i != 1; --i) {
          if (a[i] == T()) continue;
          result += detail::innerCoefficient(a[i]) + "x^" +
            boost::lexical_cast<std::string>(i);
        }

        // special case: just "x"
        if (a[i] != T()) {
          result += detail::innerCoefficient(a[i]) + "x";
        }
      }
      --i;

      // last should be just the constant
      if (a[i] != T()) {
        result += detail::innerCoefficient(a[i],detail::ShowOnes);
      }
    }
    
    return result;
  } // string formula_format(polynomial<T> const &a)
}

#endif
