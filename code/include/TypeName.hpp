/**
 * Copyright (C) 2017-2018
 * Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 *
 * @author Pablo Alvarado
 * @date   18.08.2018
 */

#ifndef ANPI_TYPE_NAME_HPP
#define ANPI_TYPE_NAME_HPP

#include <string>
#include <typeinfo>

namespace anpi {

  namespace detail {
    /// Demangle the type name 
    std::string demangle(const char* name);
  }

  /// Return the name of the given type
  template <class T>
  std::string typeName() {
    static T t;
    return detail::demangle(typeid(t).name());
  }
  
}

#endif
