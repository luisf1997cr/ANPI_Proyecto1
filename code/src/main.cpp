/**
 * Copyright (C) 2018
 * Área Académica de Ingeniería en Computadoras, TEC, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 *
 * @Author: Pablo Alvarado
 * @Date  : 24.03.2018
 */

#include <cstdlib>
#include <string>
#include <iostream>
#include <algorithm>

#include <boost/program_options.hpp>
#include <boost/type_traits/is_complex.hpp>

#include <PolynomialParser.hpp>
#include <PolynomialFormulaFormat.hpp>
#include <Exception.hpp>
#include <TypeName.hpp>

#include <Muller.hpp>
#include <JenkinsTraub.hpp>

/// Allowed types for roots and coefficients
enum TypesEnum {
  Float,
  Double,
  FComplex,
  DComplex
};

/// Allowed root finding methods
enum MethodsEnum {
  Muller,
  JenkinsTraub
};


/**
 * Pack all configuration options into one struct
 */
struct Config {
  /// Initialize with default values
  Config()
    : method(Muller)
    , polish(anpi::DoNotPolish)
    , start(0.0,0.0) {}
  
  /// Root finder method
  MethodsEnum method;
  /// Polishing roots policy
  anpi::PolishEnum polish;
  /// Starting value
  std::complex<double> start;
};

/// Convert the given string to lowercase
std::string tolower(std::string s) {
  std::transform(s.begin(),s.end(),s.begin(),
                 [](unsigned char c) -> unsigned char {
                   return std::tolower(c);
                 });
  return s;
}


/// Convert a string into the allowed type
TypesEnum as_type(const std::string& str) {
  std::string lower(tolower(str));

  if (str.find("double") != std::string::npos) {
    return Double;
  }
  if (str.find("float") != std::string::npos) {
    return Float;
  }
  if (str.find("fcomplex") != std::string::npos) {
    return FComplex;
  }
  return DComplex; // if unknown, use the most comprehensive one!
}

/// Cast complex to another complex or to a floating point type
template<class T>
inline typename std::enable_if<std::is_floating_point<T>::value,T>::type
castComplex(const std::complex<double>& a) {
  if (a.imag() == 0.0) {
    return static_cast<T>(a.real());
  }
  throw anpi::Exception("Cannot cast imaginary number into a real one");
}

template<class T>
typename std::enable_if<boost::is_complex<T>::value,T>::type
castComplex(const std::complex<double>& a) {
  typedef typename T::value_type type;
  return T(static_cast<type>(a.real()),
           static_cast<type>(a.imag()));
}


namespace bmt=boost::math::tools; // for polynomial

/**
 * Simple method to show all roots
 */
template<class T>
void showRoots(const std::vector<T>& roots) {
  std::cout << roots.size() << " roots of type "
            << anpi::typeName<T>() << " found:\n";
  for (auto& it : roots) {
    std::cout << "  " << ( anpi::detail::is_real(it) ? std::real(it) : it )
              << std::endl;
  }
}


/**
 * Call the proper method with the proper options
 */
template<class CT,class RT>
void findRoots(const Config& config,
               const std::string& polyStr) {

  // Store the roots in here
  std::vector<RT> roots;
  
  // First, convert the string to a polynomial
  bmt::polynomial<CT> poly = anpi::parsePolynomial<CT>(polyStr);

  // Show the polynomial just for reference
  std::cout << "Finding roots of type "
            << anpi::typeName<RT>() << " for:\n  "
            << anpi::polynomialFormulaFormat(poly)
            << std::endl;

  switch(config.method) {
  case Muller:
    std::cout << "with the Muller method." << std::endl;
    anpi::muller(poly,roots,config.polish,castComplex<RT>(config.start));
    break;
  case JenkinsTraub:
    std::cout << "with the Jenkins-Traub method." << std::endl;
    anpi::jenkinsTraub(poly,roots);
    break;
  default:
    throw anpi::Exception("Unknown method selected");
  }

  std::cout << "Type of coefficients: " << anpi::typeName<CT>() << std::endl;

  if (roots.empty()) {
    throw anpi::Exception("No roots found.");
  }

  showRoots(roots);
}

template<class CT>
void dispatchRootType(const Config& config,
                      const std::string& poly,
                      const TypesEnum rootType) {

  // Awful dispatch... there's no other (simpler) way to combine runtime
  // and compile-time decisions
  switch(rootType) {
  case Float: {
    findRoots<CT,float>(config,poly);
  } break;
  case Double: {
    findRoots<CT,double>(config,poly);
  } break;
  case FComplex: {
    findRoots< CT,std::complex<float> >(config,poly);
  } break;
  case DComplex: {
    findRoots< CT,std::complex<double> >(config,poly);
  } break;
  default:
    throw anpi::Exception("Unexpected root type");
  }
}

/**
 * Dispatch to the proper instantiation of the methods
 */
void dispatch(const Config& config,
              const std::string& poly,
              const TypesEnum coefType,
              const TypesEnum rootType) {

  switch(coefType) {
  case Float:
    dispatchRootType<float>(config,poly,rootType);
    break;
  case Double:
    dispatchRootType<double>(config,poly,rootType);
    break;
  case FComplex:
    dispatchRootType< std::complex<float> >(config,poly,rootType);
    break;
  case DComplex:
    dispatchRootType< std::complex<double> >(config,poly,rootType);
    break;
  default:
    throw anpi::Exception("Unexpected coefficient type");
  }
}

namespace po = boost::program_options;

////////////////////////////////////////////////////////////////////////////
//  Main program
////////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[]) {

  try {
    po::options_description desc("Allowed options");
    desc.add_options()
      ("poly,e",po::value<std::string>()->default_value(""),
       "polynomial expression")
      ("coefficients,c",po::value<std::string>()->default_value("double"),
       "force the given type for the coefficients")
      ("roots,r",po::value<std::string>()->default_value("double"),
       "expect the given type for the roots")
      ("muller,m","use the Muller method to find the roots (default)")
      ("jenkinstraub,j","use the Jenkins-Traub method to find the roots")
      ("polish,p","polish the roots")
      ("help,h", "produce help message")
      ;

    po::variables_map vm;

    po::store(parse_command_line(argc, argv, desc), vm);
    
    if (vm.count("help")) {
      std::cout
        << desc << "\n\n";
      std::cout
        << "For the coefficient and root types, use one of the following:\n"
           "  float    real numbers with single precision \n"
           "  double   real numbers with double precision \n"
           "  fcomplex complex numbers with float components\n"
           "  dcomplex complex numbers with double components\n\n"
        << "The polynomial formulae are composed by one or more \n"
           "polynomial terms:\n"
           "  polynomial  := term [ {'+'|'-'} term ]*\n"
           "  term        := coefficient 'x^' exponent\n"
           "  coefficient := double | complex\n"
           "  complex     := { '(' double ',' double ')' } | (double)\n"
           "  exponent    := unsigned integral\n\n"
        << "Examples of valid polynomials:\n"
           "  x^3 + 2x + 1\n"
           "  -5x^4 + 2.5x^2 + x\n"
           "  -3x + 5x^4 + 1 -2x^2\n"
           "  (0,1)x^4 + (5,2)x^2 + 1.5x + (1.5,2)\n"
           "The last example has some complex coefficients"
        << std::endl;
      
      return EXIT_SUCCESS;
    }

    po::notify(vm);

    std::string poly = vm["poly"].as<std::string>();

    // Default values
    TypesEnum coefType(Double),rootType(Double);
    Config config;
    config.method = Muller;
    config.start  = std::complex<double>(0.0,0.0);
    
    if (vm.count("coefficients")) {
      coefType = as_type(vm["coefficients"].as<std::string>());
    }
    
    if (vm.count("roots")) {
      rootType = as_type(vm["roots"].as<std::string>());
    }

    if (vm.count("muller")) {
      config.method = Muller;
    }

    if (vm.count("jenkinstraub")) {
      config.method = JenkinsTraub;
    }

    if (vm.count("polish")) {
      config.polish = anpi::PolishRoots;
    }
    
    // Dispatch with the proper types to call the real workers
    dispatch(config,poly,coefType,rootType);
    
  } catch(po::error& e) { 
    std::cerr << "Error:\n  " << e.what() << std::endl << std::endl; 
    return EXIT_FAILURE;
  } catch(std::exception& e)  { 
    std::cerr << "Unhandled Exception reached the top of main:\n  " 
              << e.what() << "\nApplication will now exit" << std::endl; 
    return EXIT_FAILURE;
  } 
  
  return EXIT_SUCCESS;
}
