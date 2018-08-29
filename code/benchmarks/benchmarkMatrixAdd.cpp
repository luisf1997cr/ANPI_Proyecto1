/**
 * Copyright (C) 2017-2018
 * Área Académica de Ingeniería en Computadoras, TEC, Costa Rica
 *
 * This file is part of the CE3102 Numerical Analysis lecture at TEC
 *
 * @author Pablo Alvarado
 * @date   29.12.2017
 */


#include <boost/test/unit_test.hpp>


#include <iostream>
#include <exception>
#include <cstdlib>
#include <complex>

/**
 * Unit tests for the matrix class
 */
#include "benchmarkFramework.hpp"
#include "Matrix.hpp"
#include "Allocator.hpp"

BOOST_AUTO_TEST_SUITE( Matrix )

/// Benchmark for addition operations
template<typename T>
class benchAdd {
protected:
  /// Maximum allowed size for the square matrices
  const size_t _maxSize;

  /// A large matrix holding 
  anpi::Matrix<T> _data;

  /// State of the benchmarked evaluation
  anpi::Matrix<T> _a;
  anpi::Matrix<T> _b;
  anpi::Matrix<T> _c;
public:
  /// Construct
  benchAdd(const size_t maxSize)
    : _maxSize(maxSize),_data(maxSize,maxSize,anpi::DoNotInitialize) {

    size_t idx=0;
    for (size_t r=0;r<_maxSize;++r) {
      for (size_t c=0;c<_maxSize;++c) {
        _data(r,c)=idx++;
      }
    }
  }

  /// Prepare the evaluation of given size
  void prepare(const size_t size) {
    assert (size<=this->_maxSize);
    this->_a=std::move(anpi::Matrix<T>(size,size,_data.data()));
    this->_b=this->_a;
  }
};

/// Provide the evaluation method for in-place addition 
template<typename T>
class benchAddInPlaceFallback : public benchAdd<T> {
public:
  /// Constructor
  benchAddInPlaceFallback(const size_t n) : benchAdd<T>(n) { }
  
  // Evaluate add in-place
  inline void eval() {
    anpi::fallback::add(this->_a,this->_b);
  }
};

/// Provide the evaluation method for on-copy addition 
template<typename T>
class benchAddOnCopyFallback : public benchAdd<T> {
public:
  /// Constructor
  benchAddOnCopyFallback(const size_t n) : benchAdd<T>(n) { }
  
  // Evaluate add on-copy
  inline void eval() {
    anpi::fallback::add(this->_a,this->_b,this->_c);
  }
};

/// Provide the evaluation method for in-place addition 
template<typename T>
class benchAddInPlaceSIMD : public benchAdd<T> {
public:
  /// Constructor
  benchAddInPlaceSIMD(const size_t n) : benchAdd<T>(n) { }
  
  // Evaluate add in-place
  inline void eval() {
    anpi::simd::add(this->_a,this->_b);
  }
};

/// Provide the evaluation method for on-copy addition 
template<typename T>
class benchAddOnCopySIMD : public benchAdd<T> {
public:
  /// Constructor
  benchAddOnCopySIMD(const size_t n) : benchAdd<T>(n) { }
  
  // Evaluate add on-copy
  inline void eval() {
    anpi::simd::add(this->_a,this->_b,this->_c);
  }
};

/**
 * Instantiate and test the methods of the Matrix class
 */
BOOST_AUTO_TEST_CASE( Add ) {

  std::vector<size_t> sizes = {  24,  32,  48,  64,
                                 96, 128, 192, 256,
                                384, 512, 768,1024,
                               1536,2048,3072,4096};

  const size_t n=sizes.back();
  const size_t repetitions=100;
  std::vector<anpi::benchmark::measurement> times;

  {
    benchAddOnCopyFallback<float>  baoc(n);

    // Measure on-copy add
    ANPI_BENCHMARK(sizes,repetitions,times,baoc);
    
    ::anpi::benchmark::write("add_on_copy_float_fb.txt",times);
    ::anpi::benchmark::plotRange(times,"On-copy (float) fallback","r");
  }

  {
    benchAddOnCopySIMD<float>  baoc(n);

    // Measure on-copy add
    ANPI_BENCHMARK(sizes,repetitions,times,baoc);
    
    ::anpi::benchmark::write("add_on_copy_float_simd.txt",times);
    ::anpi::benchmark::plotRange(times,"On-copy (float) simd","g");
  }
  
  {
    benchAddInPlaceFallback<float> baip(n);

    // Measure in place add
    ANPI_BENCHMARK(sizes,repetitions,times,baip);

    ::anpi::benchmark::write("add_in_place_float_fb.txt",times);
    ::anpi::benchmark::plotRange(times,"In-place (float) fallback","b");
  }

  {
    benchAddInPlaceSIMD<float> baip(n);

    // Measure in place add
    ANPI_BENCHMARK(sizes,repetitions,times,baip);

    ::anpi::benchmark::write("add_in_place_float_simd.txt",times);
    ::anpi::benchmark::plotRange(times,"In-place (float) simd","m");
  }
  
#if 0
  
  {
    benchAddOnCopy<double>  baoc(n);

    // Measure on-copy add
    ANPI_BENCHMARK(sizes,repetitions,times,baoc);
    
    ::anpi::benchmark::write("add_on_copy_double.txt",times);
    ::anpi::benchmark::plotRange(times,"On-copy (double)","g");
  }

  {
    benchAddInPlace<double> baip(n);

    // Measure in place add
    ANPI_BENCHMARK(sizes,repetitions,times,baip);

    ::anpi::benchmark::write("add_in_place_double.txt",times);
    ::anpi::benchmark::plotRange(times,"In-place (double)","m");
  }
#endif
  
  ::anpi::benchmark::show();
}
  
BOOST_AUTO_TEST_SUITE_END()
