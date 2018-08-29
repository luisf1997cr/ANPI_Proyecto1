/**
 * Copyright (C) 2017-2018
 * Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 *
 * @author Pablo Alvarado
 * @date   15.12.2017
 */

#ifndef ANPI_HAS_TYPE_HPP
#define ANPI_HAS_TYPE_HPP

/**
 * Generates a metafunction to check if a class has an internal type.
 *
 * For example:
 *
 * \code
 * GENERATE_HAS_TYPE(foobar)
 * \endcode
 *
 * generates a metafunction, such that \c has_type_foobar<T>::value
 * is true iff T::foobar is an existing type, or false otherwise.
 */
#define GENERATE_HAS_TYPE(Type)                                         \
  namespace detail {                                                    \
    template < class T >                                                \
    class HasType_##Type {                                              \
    private:                                                            \
      /* Check based on these types' sizes */                           \
      using Yes = char[2];                                              \
      using No  = char[1];                                              \
                                                                        \
      /* Fallback defines the type to cause ambiguity later */          \
      struct Fallback { struct Type { }; };                             \
      /* Derived class will have twice the Type except if NOT defined*/ \
      struct Derived : T, Fallback { };                                 \
                                                                        \
      /* SFINAE hinders instantiation of this if ambiguous */           \
      template < class U >                                              \
      static No& test ( typename U::Type* );                            \
      /* SFINAE falls back to this in case of ambiguity */              \
      template < typename U >                                           \
      static Yes& test ( U* );                                          \
    public:                                                             \
      static constexpr bool value =                                     \
        sizeof(test<Derived>(nullptr)) == sizeof(Yes);                  \
    };                                                                  \
  } /* namespace detail */                                              \
                                                                        \
  template < class T >                                                  \
  struct has_type_##Type                                                \
    : public std::integral_constant<bool,                               \
                                    detail::HasType_##Type<T>::value> { \
  }
#endif
