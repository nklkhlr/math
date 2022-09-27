//  (C) Copyright Nick Thompson 2018.
//  (C) Copyright Matt Borland 2020.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_MULTIVARIATE_STATISTICS_HPP
#define BOOST_MATH_MULTIVARIATE_STATISTICS_HPP

#include <boost/math/statistics/detail/multivariate_single_pass.hpp>

namespace boost::math::statistics {

template<class ExecutionPolicy, class ForwardIterator>
inline auto covariance(ExecutionPolicy &&exec, ForwardIterator first, ForwardIterator last)
{
    using Real = typename std::iterator_traits<ForwardIterator>::value_type;

    if constexpr (std::is_integral_v<Real>)
    {
        if constexpr (std::is_same_v<std::remove_reference_t<decltype(exec)>, decltype(std::execution::seq)>)
        {
           return detail::variance_sequential_impl<std::vector<std::vector<double>>>(first, last));
        }
        else
        {
            return detail::variance_parallel_impl<std::vector<std::vector<double>>>(first, last));
        }
    }
    else
    {
        if constexpr (std::is_same_v<std::remove_reference_t<decltype(exec)>, decltype(std::execution::seq)>)
        {
           return detail::variance_sequential_impl<std::vector<std::vector<Real>>>(first, last));
        }
        else
        {
            return detail::variance_parallel_impl<std::vector<std::vector<Real>>>(first, last));
        }

    }
}

template<class ExecutionPolicy, class Container>
inline auto covariance(ExecutionPolicy &&exec, Container const &v) {
    return covariance(exec, std::cbegin(v), std::cend(v));
}

template<class ForwardIterator>
inline auto covariance(ForwardIterator first, ForwardIterator last) {
    return covariance(std::execution::seq, first, last);
}

template<class Container>
inline auto covariance(Container const &v) {
    return covariance(std::execution::seq, std::cbegin(v), std::cend(v));
}

}

#endif //BOOST_MATH_MULTIVARIATE_STATISTICS_HPP
