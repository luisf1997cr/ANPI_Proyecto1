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

#include <boost/type_traits/is_complex.hpp>
#include <boost/math/tools/polynomial.hpp>


namespace anpi {
  
  /// Enumerator makes explicit polishing roots
  enum PolishEnum {
    DoNotPolish,
    PolishRoots
  };

  namespace bmt=boost::math::tools; // for polynomial
  
  /**
   * Compute the roots of the given polynomial using the Muller method.
   * @param[in] poly polynomial to be analyzed for roots
   * @param[out] roots all roots found
   * @param[in] start initial point for finding the roots
   * @param[in] polish indicate if polishing is needed or not.
   *
   * @return the number of roots found
   */
  template<class T,class U>
  void muller(const bmt::polynomial<T>& poly,
              std::vector<U>& roots,
              const PolishEnum polish=DoNotPolish,
              const U start=U()) {
    
    static_assert(std::is_floating_point<T>::value ||
                  boost::is_complex<T>::value,
                  "T must be floating point or complex");
   static_assert(std::is_floating_point<U>::value ||
                  boost::is_complex<U>::value,
                  "U must be floating point or complex");

//maximun allowed iterations
const unsigned int MAX_ITERATIONS= 10000;
const T eps= std::numeric_limits<T>::epsilon();

//simple Muller

//choosing 3 equally distant points to start approximation
//the choice to make the spacing 1 is arbitrary
U x1=start, x2= start+1,x3= start+2 ;

//variables to hold calculations
U q, A, B, C, xi;

//variables to hold the evaluated polynomials
U f1,f2,f3;

size_t psize= bmt::poly.size();
for(int i; i<MAX_ITERATIONS;i++){

  //evaluate polynomials at the three given approximations 
f1= bmt::evaluate_polynomial<T,U>(poly.data(),x1, psize);
f2= bmt::evaluate_polynomial<T,U>(poly.data(),x2, psize);
f3= bmt::evaluate_polynomial<T,U>(poly.data(),x3, psize);

//after evaluating the functions we check if we have arrived at a 0
if(abs(f1)<eps || abs(f1)-abs(f2)<eps){//*****************revisar, creo q esta incorrecto
return x1;

}

//calculating Mullers formulas
q=(x1-x2)/(x1-x2);
A=q*f1-q*(1+q)*f2+q*q*f3;
B=(2*q+1)*f1-(1+q)*(1+q)*f2+q*q*f3;
C=(1+q)*f1;

//calculate new approximation
xi= B>0 ? x1-(x1-x2)*(2*C)/(B + std::sqrt(B*B-4*A*C)) : x1-(x1-x2)*(2*C)/(B - std::sqrt(B*B-4*A*C));

//set values for next iteration
x3=x2;
x2=x1;
x1=xi;


}//end for
//return 

    throw Exception("No root found");
  }
}


#endif
