/*
 * Copyright (C) 2017-2018
 * Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 *
 * @Author: Pablo Alvarado
 * @Date:   28.12.2017
 */

#ifndef ANPI_MATRIX_ARITHMETIC_HPP
#define ANPI_MATRIX_ARITHMETIC_HPP

#include "Intrinsics.hpp"
#include <type_traits>

namespace anpi
{
  namespace fallback {
    /*
     * Sum
     */

    // Fallback implementation
    
    // In-copy implementation c=a+b
    template<typename T,class Alloc>
    inline void add(const Matrix<T,Alloc>& a,
                    const Matrix<T,Alloc>& b,
                    Matrix<T,Alloc>& c) {

      assert( (a.rows() == b.rows()) &&
              (a.cols() == b.cols()) );

      const size_t tentries = a.rows()*a.dcols();
      c.allocate(a.rows(),a.cols());
      
      T* here        = c.data();
      T *const end   = here + tentries;
      const T* aptr = a.data();
      const T* bptr = b.data();

      for (;here!=end;) {
        *here++ = *aptr++ + *bptr++;
      }
    }

    // In-place implementation a = a+b
    template<typename T,class Alloc>
    inline void add(Matrix<T,Alloc>& a,
                    const Matrix<T,Alloc>& b) {

      assert( (a.rows() == b.rows()) &&
              (a.cols() == b.cols()) );

      const size_t tentries = a.rows()*a.dcols();
      
      T* here        = a.data();
      T *const end   = here + tentries;
      
      const T* bptr = b.data();

      for (;here!=end;) {
        *here++ += *bptr++;
      }
    }


    /*
     * Subtraction
     */

    // Fall back implementations

    // In-copy implementation c=a-b
    template<typename T,class Alloc>
    inline void subtract(const Matrix<T,Alloc>& a,
                         const Matrix<T,Alloc>& b,
                         Matrix<T,Alloc>& c) {

      assert( (a.rows() == b.rows()) &&
              (a.cols() == b.cols()) );

      const size_t tentries = a.rows()*a.dcols();
      c.allocate(a.rows(),a.cols());
      
      T* here        = c.data();
      T *const end   = here + tentries;
      const T* aptr = a.data();
      const T* bptr = b.data();

      for (;here!=end;) {
        *here++ = *aptr++ - *bptr++;
      }
    }

    // In-place implementation a = a-b
    template<typename T,class Alloc>
    inline void subtract(Matrix<T,Alloc>& a,
                         const Matrix<T,Alloc>& b) {

      assert( (a.rows() == b.rows()) &&
              (a.cols() == b.cols()) );
      
      const size_t tentries = a.rows()*a.dcols();
      
      T* here        = a.data();
      T *const end   = here + tentries;
      
      const T* bptr = b.data();

      for (;here!=end;) {
        *here++ -= *bptr++;
      }
    }

  } // namespace fallback


  namespace simd
  {
    /*
     * Sum
     */

    /*
     * The following code exemplifies how to manually accelerate code using
     * SIMD instructions.  However, for the simple element-wise algorithms
     * like sum or subtraction, modern compilers can automatically vectorize
     * the code, as the benchmarks show.
     */


    /// We wrap the intrinsics methods to be polymorphic versions
    template<typename T,class regType>
    regType mm_add(regType,regType); // We don't implement this to cause, at
                                     // least, a linker error if this version is
                                     // used.
    //{
    // Generic function should never be called.
    // If it is called, then some SIMD chaos is going on...
    
    // A way to cause a compile time error would be better
    // throw std::bad_function_call();
    // return regType();
    //}
    
#ifdef __AVX512F__
    template<>
    inline __m512d __attribute__((__always_inline__))
    mm_add<double>(__m512d a,__m512d b) {
      return _mm512_add_pd(a,b);
    }
    template<>
    inline __m512 __attribute__((__always_inline__))
    mm_add<float>(__m512 a,__m512 b) {
      return _mm512_add_ps(a,b);
    }
    template<>
    inline __m512i __attribute__((__always_inline__))
    mm_add<uint64_t>(__m512i a,__m512i b) {
      return _mm512_add_epi64(a,b);
    }
    template<>
    inline __m512i __attribute__((__always_inline__))
    mm_add<int64_t>(__m512i a,__m512i b) {
      return _mm512_add_epi64(a,b);
    }
    template<>
    inline __m512i __attribute__((__always_inline__))
    mm_add<uint32_t>(__m512i a,__m512i b) {
      return _mm512_add_epi32(a,b);
    }
    template<>
    inline __m512i __attribute__((__always_inline__))
    mm_add<int32_t>(__m512i a,__m512i b) {
      return _mm512_add_epi32(a,b);
    }
    template<>
    inline __m512i __attribute__((__always_inline__))
    mm_add<uint16_t>(__m512i a,__m512i b) {
      return _mm512_add_epi16(a,b);
    }
    template<>
    inline __m512i __attribute__((__always_inline__))
    mm_add<int16_t>(__m512i a,__m512i b) {
      return _mm512_add_epi16(a,b);
    }
    template<>
    inline __m512i __attribute__((__always_inline__))
    mm_add<uint8_t>(__m512i a,__m512i b) {
      return _mm512_add_epi8(a,b);
    }
    template<>
    inline __m512i __attribute__((__always_inline__))
    mm_add<int8_t>(__m512i a,__m512i b) {
      return _mm512_add_epi8(a,b);
    }
#elif defined __AVX__
    template<>
    inline __m256d __attribute__((__always_inline__))
    mm_add<double>(__m256d a,__m256d b) {
      return _mm256_add_pd(a,b);
    }
    template<>
    inline __m256 __attribute__((__always_inline__))
    mm_add<float>(__m256 a,__m256 b) {
      return _mm256_add_ps(a,b);
    }
    template<>
    inline __m256i __attribute__((__always_inline__))
    mm_add<uint64_t>(__m256i a,__m256i b) {
      return _mm256_add_epi64(a,b);
    }
    template<>
    inline __m256i __attribute__((__always_inline__))
    mm_add<int64_t>(__m256i a,__m256i b) {
      return _mm256_add_epi64(a,b);
    }
    template<>
    inline __m256i __attribute__((__always_inline__))
    mm_add<uint32_t>(__m256i a,__m256i b) {
      return _mm256_add_epi32(a,b);
    }
    template<>
    inline __m256i __attribute__((__always_inline__))
    mm_add<int32_t>(__m256i a,__m256i b) {
      return _mm256_add_epi32(a,b);
    }
    template<>
    inline __m256i __attribute__((__always_inline__))
    mm_add<uint16_t>(__m256i a,__m256i b) {
      return _mm256_add_epi16(a,b);
    }
    template<>
    inline __m256i __attribute__((__always_inline__))
    mm_add<int16_t>(__m256i a,__m256i b) {
      return _mm256_add_epi16(a,b);
    }
    template<>
    inline __m256i __attribute__((__always_inline__))
    mm_add<uint8_t>(__m256i a,__m256i b) {
      return _mm256_add_epi8(a,b);
    }
    template<>
    inline __m256i __attribute__((__always_inline__))
    mm_add<int8_t>(__m256i a,__m256i b) {
      return _mm256_add_epi8(a,b);
    }
#elif  defined __SSE2__
    template<>
    inline __m128d __attribute__((__always_inline__))
    mm_add<double>(__m128d a,__m128d b) {
      return _mm_add_pd(a,b);
    }
    template<>
    inline __m128 __attribute__((__always_inline__))
    mm_add<float>(__m128 a,__m128 b) {
      return _mm_add_ps(a,b);
    }
    template<>
    inline __m128i __attribute__((__always_inline__))
    mm_add<std::uint64_t>(__m128i a,__m128i b) {
      return _mm_add_epi64(a,b);
    }
    template<>
    inline __m128i __attribute__((__always_inline__))
    mm_add<std::int64_t>(__m128i a,__m128i b) {
      return _mm_add_epi64(a,b);
    }
    template<>
    inline __m128i __attribute__((__always_inline__))
    mm_add<std::uint32_t>(__m128i a,__m128i b) {
      return _mm_add_epi32(a,b);
    }
    template<>
    inline __m128i __attribute__((__always_inline__))
    mm_add<std::int32_t>(__m128i a,__m128i b) {
      return _mm_add_epi16(a,b);
    }
    template<>
    inline __m128i __attribute__((__always_inline__))
    mm_add<std::uint16_t>(__m128i a,__m128i b) {
      return _mm_add_epi16(a,b);
    }
    template<>
    inline __m128i __attribute__((__always_inline__))
    mm_add<std::int16_t>(__m128i a,__m128i b) {
      return _mm_add_epi32(a,b);
    }
    template<>
    inline __m128i __attribute__((__always_inline__))
    mm_add<std::uint8_t>(__m128i a,__m128i b) {
      return _mm_add_epi16(a,b);
    }
    template<>
    inline __m128i __attribute__((__always_inline__))
    mm_add<std::int8_t>(__m128i a,__m128i b) {
      return _mm_add_epi32(a,b);
    }
#endif
    
    // On-copy implementation c=a+b
    template<typename T,class Alloc,typename regType>
    inline void addSIMD(const Matrix<T,Alloc>& a, 
                        const Matrix<T,Alloc>& b,
                        Matrix<T,Alloc>& c) {

      // This method is instantiated with unaligned allocators.  We
      // allow the instantiation although externally this is never
      // called unaligned
      static_assert(!extract_alignment<Alloc>::aligned ||
		    (extract_alignment<Alloc>::value >= sizeof(regType)),
		    "Insufficient alignment for the registers used");
      
      const size_t tentries = a.rows()*a.dcols();
      c.allocate(a.rows(),a.cols());

      regType* here        = reinterpret_cast<regType*>(c.data());
      const size_t  blocks = ( tentries*sizeof(T) + (sizeof(regType)-1) )/
        sizeof(regType);
      regType *const end   = here + blocks;
      const regType* aptr  = reinterpret_cast<const regType*>(a.data());
      const regType* bptr  = reinterpret_cast<const regType*>(b.data());
      
      for (;here!=end;) {
        *here++ = mm_add<T>(*aptr++,*bptr++);
      }
      
    }
       
    // On-copy implementation c=a+b for SIMD-capable types
    template<typename T,
	     class Alloc,
	     typename std::enable_if<is_simd_type<T>::value,int>::type=0>
    inline void add(const Matrix<T,Alloc>& a,
                    const Matrix<T,Alloc>& b,
                    Matrix<T,Alloc>& c) {

      assert( (a.rows() == b.rows()) &&
              (a.cols() == b.cols()) );


      if (is_aligned_alloc<Alloc>::value) {        
#ifdef __AVX512F__
        addSIMD<T,Alloc,typename avx512_traits<T>::reg_type>(a,b,c);
#elif  __AVX__
        addSIMD<T,Alloc,typename avx_traits<T>::reg_type>(a,b,c);
#elif  __SSE2__
        addSIMD<T,Alloc,typename sse2_traits<T>::reg_type>(a,b,c);
#else
        ::anpi::fallback::add(a,b,c);
#endif
      } else { // allocator seems to be unaligned
        ::anpi::fallback::add(a,b,c);
      }
    }

    // Non-SIMD types such as complex
    template<typename T,
             class Alloc,
             typename std::enable_if<!is_simd_type<T>::value,int>::type = 0>
    inline void add(const Matrix<T,Alloc>& a,
                    const Matrix<T,Alloc>& b,
                    Matrix<T,Alloc>& c) {
      
      ::anpi::fallback::add(a,b,c);
    }

    // In-place implementation a = a+b
    template<typename T,class Alloc>
    inline void add(Matrix<T,Alloc>& a,
                    const Matrix<T,Alloc>& b) {

      add(a,b,a);
    }


    /*
     * Subtraction
     */

    // Fall back implementations

    // In-copy implementation c=a-b
    template<typename T,class Alloc>
    inline void subtract(const Matrix<T,Alloc>& a,
                         const Matrix<T,Alloc>& b,
                         Matrix<T,Alloc>& c) {
      ::anpi::fallback::subtract(a,b,c);
    }

    // In-place implementation a = a-b
    template<typename T,class Alloc>
    inline void subtract(Matrix<T,Alloc>& a,
                         const Matrix<T,Alloc>& b) {

      ::anpi::fallback::subtract(a,b);
    }
  } // namespace simd


  // The arithmetic implementation (aimpl) namespace
  // dispatches to the corresponding methods
#ifdef ANPI_ENABLE_SIMD
  namespace aimpl=simd;
#else
  namespace aimpl=fallback;
#endif
  
} // namespace anpi

#endif
