/**
 * Copyright (C) 2017-2018
 * Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 *
 * @author Pablo Alvarado
 * @date   15.12.2017
 */

#ifndef ANPI_ALLOCATOR_HPP
#define ANPI_ALLOCATOR_HPP

#include <boost/align/aligned_allocator.hpp>
#include "HasType.hpp"

namespace anpi {

# if defined __AVX512F__
  static const size_t DefaultAlignment = 64;
# elif defined __AVX2__
  static const size_t DefaultAlignment = 32;
# elif defined __AVX__
  static const size_t DefaultAlignment = 32;
# elif defined __SSE4_2__
  static const size_t DefaultAlignment = 16;
# elif defined __SSE4_1__
  static const size_t DefaultAlignment = 16;
# elif defined __SSSE3__
  static const size_t DefaultAlignment = 16;
# elif defined __SSE3__
  static const size_t DefaultAlignment = 16;
# elif defined __SSE2__
  static const size_t DefaultAlignment = 16;
# elif defined __SSE__
  static const size_t DefaultAlignment = 16;
# endif  
  
  /**
   * Use the boost version of aligned_allocator
   */
  template<class T, std::size_t Align=DefaultAlignment>
  class aligned_allocator : public boost::alignment::aligned_allocator<T,Align>
  {
  public:
    /// Inherit all constructors
    using boost::alignment::aligned_allocator<T,Align>::aligned_allocator;

    /// Change the stored type
    template<class U>
    struct rebind {
      typedef aligned_allocator<U, Align> other;
    };
  };

  /**
   * This is identical to the aligned_allocator, but additionally
   * it exposes the used Alignment in a static attribute
   */
  template<class T, std::size_t Align=DefaultAlignment>
  class aligned_row_allocator : public aligned_allocator<T,Align> {
  public:
    /// Inherit all constructors
    using aligned_allocator<T,Align>::aligned_allocator;

    /// Change the stored type
    template<class U>
    struct rebind {
      typedef aligned_row_allocator<U, Align> other;
    };

    /// Type to identify this as a row-aligned allocator
    typedef std::true_type row_aligned;
  };


  /**
   * Check if a class is an aligned_allocator or an aligned_row_allocator
   */
  template<class Alloc>
  struct is_aligned_alloc {
    static const bool value = false;
  };

  // Specialization for the aligned allocator
  template<typename T, std::size_t A>
  struct is_aligned_alloc< anpi::aligned_allocator<T,A> > {
    static const bool value = true;
  };

  // Specialization for the aligned_row_allocator
  template<typename T, std::size_t A>
  struct is_aligned_alloc< anpi::aligned_row_allocator<T,A> > {
    static const bool value = true;
  };

  /**
   * Create metafunction has_type_row_aligned<T>
   */
  GENERATE_HAS_TYPE(row_aligned);
  
  /**
   * Extract alignment of an allocator
   *
   * We assume that the first template type of the allocator corresponds
   * to the type of the elements of the allocated space.
   *
   * If a second integral template parameter is provided, it will be assume
   * to correspond to the alignment.
   */
  template<class Alloc>
  struct extract_alignment {
    static constexpr size_t value       = sizeof(typename Alloc::value_type);
    static constexpr bool   aligned     = false;
    static constexpr bool   row_aligned = false;
  };

  /**
   * Specialization covers both aligned allocators
   *
   * We assume that the second integral parameter of the template type
   * corresponds to the alignment parameter.
   */
  template< template<typename,size_t> class Alloc,typename T,size_t Align >
  struct extract_alignment< Alloc<T,Align> > {
    static constexpr size_t value       = Align;
    static constexpr bool   aligned     = true;
    static constexpr bool   row_aligned =
      has_type_row_aligned<Alloc<T,Align> >::value;
  };
  
}

#endif
