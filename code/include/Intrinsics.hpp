/*
 * Copyright (C) 2017-2018
 * Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 *
 * @Author: Pablo Alvarado
 * @Date:   28.12.2017
 */

#ifndef ANPI_INTRINSICS_HPP
#define ANPI_INTRINSICS_HPP

#include <cstdint>

/*
 * Include the proper intrinsics headers for the current architecture
 */

#if defined(_MSC_VER)
   /* Microsoft C/C++-compatible compiler */
#  include <intrin.h>
#elif defined(__GNUC__)
#  if (defined(__x86_64__) || defined(__i386__))
     /* GCC-compatible compiler, targeting x86/x86-64 */
#    include <x86intrin.h>
#  elif defined(__ARM_NEON__)
     /* GCC-compatible compiler, targeting ARM with NEON */
#    include <arm_neon.h>
#  elif defined(__IWMMXT__)
     /* GCC-compatible compiler, targeting ARM with WMMX */
#    include <mmintrin.h>
#  elif defined(__xlC__)) && (defined(__VEC__) || defined(__ALTIVEC__))
     /* XLC or GCC-compatible compiler, targeting PowerPC with VMX/VSX */
#    include <altivec.h>
#  elif defined(__GNUC__) && defined(__SPE__)
     /* GCC-compatible compiler, targeting PowerPC with SPE */
#    include <spe.h>
#  endif
#endif

template <typename T>
struct is_simd_type {
  static constexpr bool value =
    std::is_same<T,double>::value        ||
    std::is_same<T,float>::value         ||
    std::is_same<T,std::int64_t>::value  ||
    std::is_same<T,std::uint64_t>::value ||    
    std::is_same<T,std::int32_t>::value  ||
    std::is_same<T,std::uint32_t>::value ||
    std::is_same<T,std::int16_t>::value  ||
    std::is_same<T,std::uint16_t>::value ||
    std::is_same<T,std::int8_t>::value   ||
    std::is_same<T,std::uint8_t>::value;
};


#ifdef __AVX512F__
template<typename T> struct avx512_traits { };
template<> struct avx512_traits<double> { typedef __m512d reg_type; };
template<> struct avx512_traits<float> { typedef __m512 reg_type; };
template<> struct avx512_traits<int64_t> { typedef __m512i reg_type; };
template<> struct avx512_traits<uint64_t> { typedef __m512i reg_type; };
template<> struct avx512_traits<int32_t> { typedef __m512i reg_type; };
template<> struct avx512_traits<uint32_t> { typedef __m512i reg_type; };
template<> struct avx512_traits<int16_t> { typedef __m512i reg_type; };
template<> struct avx512_traits<uint16_t> { typedef __m512i reg_type; };
template<> struct avx512_traits<int8_t> { typedef __m512i reg_type; };
template<> struct avx512_traits<uint8_t> { typedef __m512i reg_type; };
#endif

#ifdef __AVX__
template<typename T> struct avx_traits { };
template<> struct avx_traits<double> { typedef __m256d reg_type; };
template<> struct avx_traits<float> { typedef __m256 reg_type; };
template<> struct avx_traits<int64_t> { typedef __m256i reg_type; };
template<> struct avx_traits<uint64_t> { typedef __m256i reg_type; };
template<> struct avx_traits<int32_t> { typedef __m256i reg_type; };
template<> struct avx_traits<uint32_t> { typedef __m256i reg_type; };
template<> struct avx_traits<int16_t> { typedef __m256i reg_type; };
template<> struct avx_traits<uint16_t> { typedef __m256i reg_type; };
template<> struct avx_traits<int8_t> { typedef __m256i reg_type; };
template<> struct avx_traits<uint8_t> { typedef __m256i reg_type; };
#endif

#ifdef __SSE2__
template<typename T> struct sse2_traits { };
template<> struct sse2_traits<double> { typedef __m128d reg_type; };
template<> struct sse2_traits<float> { typedef __m128 reg_type; };
template<> struct sse2_traits<int64_t> { typedef __m128i reg_type; };
template<> struct sse2_traits<uint64_t> { typedef __m128i reg_type; };
template<> struct sse2_traits<int32_t> { typedef __m128i reg_type; };
template<> struct sse2_traits<uint32_t> { typedef __m128i reg_type; };
template<> struct sse2_traits<int16_t> { typedef __m128i reg_type; };
template<> struct sse2_traits<uint16_t> { typedef __m128i reg_type; };
template<> struct sse2_traits<int8_t> { typedef __m128i reg_type; };
template<> struct sse2_traits<uint8_t> { typedef __m128i reg_type; };
#endif
  
  
#endif
