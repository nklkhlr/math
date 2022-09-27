//  (C) Copyright Nick Thompson 2019.
//  (C) Copyright Matt Borland 2021.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_HOTELLING_T2_HPP
#define BOOST_MATH_HOTELLING_T2_HPP

#include <cmath>
#include <cstddef>
#include <iterator>
#include <utility>
#include <type_traits>
#include <vector>
#include <stdexcept>
#include <boost/math/tools/assert.hpp>
#include <boost/math/statistics/detail/single_pass.hpp>
#include <boost/math/statistics/detail/multivariate_single_pass.hpp>
#include <boost/math/distributions/fisher_f.hpp>

namespace boost { namespace math { namespace statistics { namespace detail {

template<typename SizeType, typename DataType>
std::pair<DataType, std::size_t> hotelling_t2_var_equal_statistics(SizeType n_1, SizeType n_2, SizeType p, std::vector<DataType> mean_1, std::vector<DataType> mean_2, std::vector<std::vector<DataType>> cov_1, std::vector<std::vector<DataType>> cov_2)
{
	DataType test_statistic;
	std::size_t dof;
    // TODO
	return std::make_pair(test_statistic, dof);
}

template<typename SizeType, typename DataType>
 std::pair<DataType, std::size_t> hotelling_t2_var_unequal_statistics(SizeType n_1, SizeType n_2, SizeType p, std::vector<DataType> mean_1, std::vector<DataType> mean_2, std::vector<std::vector<DataType>> cov_1, std::vector<std::vector<DataType>> cov_2)
{
	DataType test_statistic;
	std::size_t dof;
    // TODO
	return std::make_pair(test_statistic, dof);
}

template<typename ReturnType, typename ForwardIterator>
ReturnType hotelling_t2_test_impl(ForwardIterator begin_1, ForwardIterator end_1, ForwardIterator begin_2, ForwardIterator end_2, bool var_equal=true) {
    using Real = typename std::tuple_element<0, ReturnType>::type;
    using std::sqrt;
	// number of samples
    auto n1 = std::distance(begin_1, end_1);
    auto n2 = std::distance(begin_2, end_2);

    // number of features (using the sample only)
    auto p = std::distance(std::begin(*begin_1), std::end(*begin_1));
    auto p2 = std::distance(std::begin(*begin_2), std::end(*begin_2));
    BOOST_MATH_ASSERT_MSG(p == p2, "Both groups must have the same number of features.");

	// computing mean vectors
    std::vector<Real> mean_1(p), mean_2(p);
    // TODO: move means to multivariate_single_pass
    std::size_t i = 0;
    while (begin_1 != end_1) {
        mean_1[i] = mean_sequential_impl(*begin_1.begin(), *begin_1.end());
        ++i;
        ++begin_1;
    }
    i = 0;
    while (begin_2 != end_2) {
        mean_2[i] = mean_sequential_impl(*begin_2.begin(), *begin_2.end());
        ++i;
        ++begin_2;
    }

    // TODO: test covariances
	std::vector<std::vector<Real>> cov_1, cov_2;
    cov_1 = covariance(begin_1, end_1);
    cov_2 = covariance(begin_2, end_2);

	// computing test statistics (t2)
	std::pair<Real, std::size_t> t2_stats;
	if (var_equal) t2_stats = hotelling_t2_var_equal_statistics(n1, n2, p, mean_1, mean_2, cov_1, cov_2);
	else t2_stats = hotelling_t2_var_unequal_statistics(n1, n2, p, mean_1, mean_2, cov_1, cov_2);
    // extracting statistics
    Real test_statistic = std::get<0>(t2_stats);
	std::size_t dof = std::get<1>(t2_stats);

    // computing p-value
    auto f_dist = boost::math::fisher_f_distribution<Real, no_promote_policy>(p, dof);
    Real pvalue = boost::math::pdf<Real>(f_dist, test_statistic);

    return std::make_pair(test_statistic, pvalue);
}

} // namespace detail

// TODO: add var_equal
template<typename ForwardIterator, typename Real = typename std::iterator_traits<ForwardIterator>::value_type,
        typename std::enable_if<std::is_integral<Real>::value, bool>::type = true>
inline auto hotelling_t2_test(ForwardIterator begin_1, ForwardIterator end_1, ForwardIterator begin_2, ForwardIterator end_2) -> std::pair<double, double>
{
    return detail::hotelling_t2_test_impl<std::pair<double, double>>(begin_1, end_1, begin_2, end_2);
}

template<typename ForwardIterator, typename Real = typename std::iterator_traits<ForwardIterator>::value_type,
        typename std::enable_if<!std::is_integral<Real>::value, bool>::type = true>
inline auto hotelling_t2_test(ForwardIterator begin_1, ForwardIterator end_1, ForwardIterator begin_2, ForwardIterator end_2) -> std::pair<Real, Real>
{
    return detail::hotelling_t2_test_impl<std::pair<Real, Real>>(begin_1, end_1, begin_2, end_2);
}

template<typename Container, typename Real = typename Container::value_type, typename std::enable_if<std::is_integral<Real>::value, bool>::type = true>
inline auto hotelling_t2_test(Container const & u, Container const & v) -> std::pair<double, double>
{
    return detail::hotelling_t2_test_impl<std::pair<double, double>>(std::begin(u), std::end(u), std::begin(v), std::end(v));
}

template<typename Container, typename Real = typename Container::value_type, typename std::enable_if<!std::is_integral<Real>::value, bool>::type = true>
inline auto hotelling_t2_test(Container const & u, Container const & v) -> std::pair<Real, Real>
{
    return detail::hotelling_t2_test_impl<std::pair<Real, Real>>(std::begin(u), std::end(u), std::begin(v), std::end(v));
}

}}} // namespace boost::math::statistics

#endif //BOOST_MATH_HOTELLING_T2_HPP
