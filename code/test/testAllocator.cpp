/**
 * Copyright (C) 2017-2018
 * Área Académica de Ingeniería en Computadoras, TEC, Costa Rica
 *
 * This file is part of the CE3102 Numerical Analysis lecture at TEC
 *
 * @author: Pablo Alvarado
 * @date:   31.12.2017
 */

#include <boost/test/unit_test.hpp>
#include <Allocator.hpp>

#define COMMA ,

BOOST_AUTO_TEST_SUITE( Allocator )

BOOST_AUTO_TEST_CASE( Allocation ) {

  {
    typedef anpi::aligned_row_allocator<float,32> alloc_type;
    alloc_type alloc;
    alloc_type::pointer ptr = alloc.allocate(1024);
    size_t ptrcst = reinterpret_cast<size_t>(ptr);
    
    BOOST_CHECK( ptrcst % 32 == 0);
    
    alloc.deallocate(ptr,1024);
  }

  {
    typedef anpi::aligned_allocator<double,64> alloc_type;
    alloc_type alloc;
    alloc_type::pointer ptr = alloc.allocate(1024);
    size_t ptrcst = reinterpret_cast<size_t>(ptr);
    
    BOOST_CHECK( ptrcst % 64 == 0);
    
    alloc.deallocate(ptr,1024);
  }
  
}

BOOST_AUTO_TEST_CASE( Checks ) {
  BOOST_CHECK(!anpi::is_aligned_alloc<std::allocator<float> >::value );
  bool val = anpi::is_aligned_alloc<anpi::aligned_allocator<float,16> >::value;
  BOOST_CHECK(val);
  val = anpi::is_aligned_alloc<anpi::aligned_row_allocator<float,16> >::value;
  BOOST_CHECK(val);

  {
    typedef anpi::extract_alignment<std::allocator<float> > extractor;
    BOOST_CHECK(extractor::value==sizeof(float) );
    BOOST_CHECK(extractor::aligned == false );
    BOOST_CHECK(extractor::row_aligned == false );
  }

  {
    typedef anpi::extract_alignment<anpi::aligned_allocator<float,16> > ext;
    BOOST_CHECK(ext::value==16 );
    BOOST_CHECK(ext::aligned == true );
    BOOST_CHECK(ext::row_aligned == false );
  }

  {
    typedef anpi::extract_alignment<anpi::aligned_row_allocator<int,32> > ext;
    BOOST_CHECK(ext::value==32);
    BOOST_CHECK(ext::aligned == true );
    BOOST_CHECK(ext::row_aligned == true );
  }

  
}

BOOST_AUTO_TEST_SUITE_END()
