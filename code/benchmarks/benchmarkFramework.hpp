/**
 * Copyright (C) 2017-2018
 * Área Académica de Ingeniería en Computadoras, TEC, Costa Rica
 *
 * This file is part of the CE3102 Numerical Analysis lecture at TEC
 *
 * @author Pablo Alvarado
 * @date   29.12.2017
 */

#ifndef ANPI_BENCHMARK_FRAMEWORK_HPP
#define ANPI_BENCHMARK_FRAMEWORK_HPP

#include <chrono>
#include <iostream>
#include <ostream>
#include <fstream>
#include <limits>

#include <Matrix.hpp>
#include <PlotPy.hpp>


namespace anpi {
  namespace benchmark {
    /**
     * Each measurement is composed by four attributes
     */
    struct measurement {
      inline measurement() : size(0u),average(0.),stddev(0.),min(0.),max(0.) {};
      
      size_t size;
      double average;
      double stddev;
      double min;
      double max;
    };

    template<typename T>
    inline T sqr(const T val) { return val*val; }
    
    /**
     * Compute measurement statistics for each size
     */
    void computeStats(const std::vector<size_t>& sizes,
                      const anpi::Matrix<std::chrono::duration<double> >& mat,
                      std::vector<measurement>& times) {

      const size_t nums = sizes.size();
      times.resize(nums);
      for ( size_t s=0;s<nums;++s ) {
        measurement& m = times[s];
        m.size = sizes[s];        

        double val = mat[s][0].count();
        m.average = val;
        m.stddev  = sqr(val);

        m.min = val;
        m.max = val;
        
        for ( size_t i=1;i<mat.cols();++i ) {
          val = (mat[s][i]-mat[s][i-1]).count();
          m.average += val;
          m.stddev  += sqr(val);
          m.min = std::min(m.min,val);
          m.max = std::max(m.max,val);
        }

        m.average /= mat.cols();
        m.stddev = std::sqrt(m.stddev/mat.cols() - sqr(m.average));
      }
    }

    /**
     * Save a file with each measurement in a row.
     *
     * The meaning of the columns is as follows:
     * # Size
     * # Average
     * # Standard deviation
     * # Minimum
     * # Maximum  
     */
    void write(std::ostream& stream,
               const std::vector<measurement>& m) {
      for (auto i : m) {
        stream << i.size    << " \t";
        stream << i.average << " \t";
        stream << i.stddev  << " \t";
        stream << i.min     << " \t";
        stream << i.max     << " \t" << std::endl;
      }
    }

    /**
     * Save a file with each measurement in a row
     */
    void write(const std::string& filename,
               const std::vector<measurement>& m) {
      std::ofstream os(filename.c_str());
      write(os,m);
      os.close();
    }

    /**
     * Plot measurements (average only)
     *
     * The meaning of the columns is as follows:
     * # Size
     * # Average
     * # Standard deviation
     * # Minimum
     * # Maximum  
     */
    void plot(const std::vector<measurement>& m,
              const std::string& legend,
              const std::string& color = "r") {
      std::vector<double> x(m.size()),y(m.size());

      for (size_t i=0;i<m.size();++i) {
        x[i]=m[i].size;
        y[i]=m[i].average;
      }

      static anpi::Plot2d<double> plotter;
      plotter.initialize(1);
      plotter.plot(x,y,legend,color);
    }

    /**
     * Plot measurements (average only)
     *
     * The meaning of the columns is as follows:
     * # Size
     * # Average
     * # Standard deviation
     * # Minimum
     * # Maximum  
     */
    void plotRange(const std::vector<measurement>& m,
                   const std::string& legend,
                   const std::string& color) {
      std::vector<double> x(m.size()),y(m.size()),miny(m.size()),maxy(m.size());

      for (size_t i=0;i<m.size();++i) {
        const measurement& mi = m[i];
        x[i]=mi.size;
        y[i]=mi.average;
        miny[i]=mi.min;
        maxy[i]=mi.max;
      }

      static anpi::Plot2d<double> plotter;
      plotter.initialize(1);
      plotter.plot(x,y,miny,maxy,legend,color);
    }
    
    void show() {
       static anpi::Plot2d<double> plotter;
       plotter.show();
    }
  } // namespace benchmark
} // namespace anpi
    
/**
 * Meassure the time for all given sizes.
 * @param sizes  vector with all sizes to be tested
 * @param rep    number of repetitions to meassure the time
 * @param times  measurement taken for each time
 *               its time must be anpi::Matrix<anpi::benchmark::measurement>
 * @param bench  benchmark instance.  See below for requirements
 *
 * The @bench is an instance of a class that must provide at least
 * the following:
 * - an void prepare(const size_t size) method, that initializes the
 *   state of the instance as required for the evaluation.  This
 *   method is called outside the performance measurements.
 * - an inline void eval() method that performs the evaluation.
 *
 * The following variables are available in PREFIX and CODE:
 * - size holds the current evaluated size (type size_t)
 * - s    holds the index of the current evaluated size
 * - i    index of the current repetition
 */
#define ANPI_BENCHMARK(sizes,rep,times,bench)                      \
{                                                                  \
  /* number of sizes to be tested */                               \
  const size_t _nums = sizes.size();                               \
                                                                   \
  typedef std::chrono::duration<double> durat;                     \
                                                                   \
  /* each row holds a particular size */                           \
  ::anpi::Matrix<durat> _mat(_nums,rep,::anpi::DoNotInitialize);   \
                                                                   \
  /* test each size */                                             \
  for (size_t s=0;s<_nums;++s ) {                                  \
    const size_t size = sizes[s];                                  \
                                                                   \
    std::cout << "Testing size " << size << std::endl;             \
                                                                   \
    durat* _row = _mat[s];                                         \
                                                                   \
    /* prefix code for initialization before measurement */        \
    bench.prepare(size);                                           \
                                                                   \
    const auto _start = std::chrono::high_resolution_clock::now(); \
                                                                   \
    for ( size_t i=0;i<rep;++i ) {                                 \
      /* the code to be benchmarked */                             \
      bench.eval();                                                \
                                                                   \
      _row[i] = std::chrono::high_resolution_clock::now()-_start;  \
    }                                                              \
  }                                                                \
                                                                   \
  ::anpi::benchmark::computeStats(sizes,_mat,times);               \
}


#endif
