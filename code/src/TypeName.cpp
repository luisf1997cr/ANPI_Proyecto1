/**
 * Copyright (C) 2017-2018
 * Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 *
 * @author Pablo Alvarado
 * @date   18.08.2018
 */

#include <string>
#include <typeinfo>

#include "TypeName.hpp"

#ifdef __GNUG__

#include <cstdlib>
#include <memory>
#include <cxxabi.h>

namespace anpi {
  namespace detail {
    /// Demangle the type name 
    std::string demangle(const char* name) {
      int status = -4; // some arbitrary value to eliminate the compiler warning

      // enable c++11 by passing the flag -std=c++11 to g++
      std::unique_ptr<char, void(*)(void*)> res {
        abi::__cxa_demangle(name, NULL, NULL, &status),
          std::free
          };
      
      return (status==0) ? res.get() : name ;
    }
  }
}

#else

namespace anpi {
  namespace detail {
    /// Does nothing if not g++
    std::string demangle(const char* name) {
      return name;
    }
  }
}

#endif


