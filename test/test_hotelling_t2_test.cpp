//  (C) Copyright Matt Borland 2021.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "math_unit_test.hpp"
#include <boost/math/statistics/hotelling_t2.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <limits>
#include <vector>
#include <random>
#include <utility>

using quad = boost::multiprecision::cpp_bin_float_quad;
using std::sqrt;

template<typename Real>
void test_hotelling_t2()
{
    std::vector<Real> m_1 {{1,2,3,4,5}, {}};
    std::vector<Real> m_2 {{2,3,4,5,6}, {}};

    std::pair<Real, Real> temp = boost::math::statistics::hotelling_t2_test(m_1, m_2);
    Real computed_statistic = std::get<0>(temp);
    Real computed_pvalue = std::get<1>(temp);
    CHECK_ULP_CLOSE(Real(1), computed_statistic, 5);
    CHECK_MOLLIFIED_CLOSE(Real(0), computed_pvalue, sqrt(std::numeric_limits<Real>::epsilon()));
}


int main()
{
    test_hotelling_t2<float>();
    test_hotelling_t2<double>();
    test_hotelling_t2<quad>();
}
