//  Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Basic sanity check that header <boost/math/special_functions/expm1.hpp>
// #includes all the files that it needs to.
//
#include <boost/math/special_functions/expm1.hpp>
//
// Note this header includes no other headers, this is
// important if this test is to be meaningful:
//
#include "test_compile_result.hpp"

void check()
{
   check_result<float>(boost::math::expm1<float>(f));
   check_result<double>(boost::math::expm1<double>(d));
   check_result<long double>(boost::math::expm1<long double>(l));
}