#pragma once
#include <algorithm>
#include <type_traits>
#include <tuple>

#include "constants.h"
#include "forward_declarations.h"
#include "gcd.h"
#include "rational.h"
#include "type_at.h"
#include "base_types.h"
#include "helpers.h"

namespace euler
{
	template<typename T, template<auto x, auto y> typename rational_type, typename ...coefficients>
	struct polynomial
	{
		// this is bad, as x and y could be anything. 
		using zero_type = polynomial<T, rational_type, typename rational_type<0, 1>::zero_type>;
		using one_type = polynomial<T, rational_type, typename rational_type<0, 1>::one_type>;

		using X = polynomial<T, rational_type, typename rational_type<0, 1>::zero_type, typename rational_type<0, 1>::one_type>;

		template<int64_t p>
		static CONST monomial() {
			if constexpr (p <= 0)
			{
				return one_type{};
			}
			else
			{
				return monomial<p - 1>() * X {};
			}
		}

		constexpr bool is_zero() const {
			return sizeof...(coefficients) == 1 && std::is_same<typename type_at<0, coefficients...>::type, typename type_at<0, coefficients...>::type::zero_type>::value;
		}

		constexpr int64_t degree() const {
			return degree_helper<rational_type, sizeof...(coefficients) - 1, coefficients...>::value();
		}

		template<int64_t index>
		CONST coeff_at() const {
			if constexpr (index < 0 || index >= sizeof...(coefficients)) {
				return typename type_at<0, coefficients...>::type::zero_type{};
			}
			else {
				return typename type_at<index, coefficients...>::type{};
			}
		}

		constexpr T operator()(const T& x) const {
			if constexpr (sizeof...(coefficients) == 0)
			{
				return constants<T>::zero;
			}
			else
			{
				return poly_eval_helper<sizeof...(coefficients) - 1, sizeof...(coefficients)>::template apply<T, coefficients...>(constants<T>::zero, x);
			}
		}

		template<typename ...Bs>
		CONST operator+(const polynomial<T, rational_type, Bs...>& b) const
		{
			return simplify(add(*this, b));
		}

		template<typename ...Bs>
		CONST operator-(const polynomial<T, rational_type, Bs...>& b) const
		{
			return simplify(sub(*this, b));
		}

		template<typename ...Bs>
		CONST operator*(const polynomial<T, rational_type, Bs...> b) const {
			return mul(*this, b);
		}

		template<auto x, auto y>
		CONST scale(const rational_type<x, y>) const {
			return mul(*this, polynomial<T, rational_type, rational_type<x, y>>{});
		}

		template<typename ...Bs>
		CONST operator/(const polynomial<T, rational_type, Bs...> b) const {
			return div(*this, b);
		}
	};

	template<typename T, template<auto, auto> typename rational_type, typename ...Ns, typename ...Ds>
	CONST div(const polynomial<T, rational_type, Ns...> n, const polynomial<T, rational_type, Ds...> d) {
		CONST q = typename polynomial<T, rational_type, typename type_at<0, Ns...>::type>::zero_type{};
		return div_helper<T, rational_type, Ns...>::template inner<Ds...>::apply(n, d, q, n);
	}

	template<typename T, template<auto x, auto y> typename rational_type, typename ...As, typename ...Bs>
	CONST add(polynomial<T, rational_type, As...> a, polynomial<T, rational_type, Bs...> b) {
		return add_low(a, b, std::make_index_sequence<std::max(sizeof...(As), sizeof...(Bs))>());
	}

	template<typename T, template<auto x, auto y> typename rational_type, typename ...As, typename ...Bs>
	CONST sub(polynomial<T, rational_type, As...> a, polynomial<T, rational_type, Bs...> b) {
		return sub_low(a, b, std::make_index_sequence<std::max(sizeof...(As), sizeof...(Bs))>());
	}

	template<typename T, template<auto x, auto y> typename rational_type, typename ...As, typename ...Bs>
	CONST mul(polynomial<T, rational_type, As...> a, polynomial<T, rational_type, Bs...> b) {
		return mul_low(a, b, std::make_index_sequence<sizeof...(As) + sizeof...(Bs) - 1>());
	}

	template<template<auto index> typename, typename, class>
	struct make_taylor_impl;

	template<template<auto index> typename coeff_at, typename T, int64_t... Is>
	struct make_taylor_impl<coeff_at, T, std::integer_sequence<int64_t, Is...>> {
		using type = polynomial<T, rational64, typename coeff_at<Is>::type...>;
	};

	// generic taylor serie, depending on coefficients
	// TODO: use that for all series below
	template<typename T, template<auto index> typename coeff_at, int64_t deg>
	using taylor = typename make_taylor_impl<coeff_at, T, std::make_integer_sequence<int64_t, deg>>::type;

	template<int64_t i>
	struct lnp1_coeff { using type = rational64<integers::alternate<i + 1>::val, i>; };

	template<>
	struct lnp1_coeff<0> { using type = rational64<0, 1>; };

	template<int64_t i, typename E = void>
	struct sin_coeff_helper {};

	template<int64_t i>
	struct sin_coeff_helper<i, typename std::enable_if<(i & 1) == 0>::type> {
		using type = rational64<0, 1>;
	};

	template<int64_t i>
	struct sin_coeff_helper<i, typename std::enable_if<(i & 1) == 1>::type> {
		using type = rational64<integers::alternate<i / 2>::val, integers::factorial<i>::val>;
	};

	template<int64_t i>
	struct sin_coeff {
		using type = typename sin_coeff_helper<i>::type;
	};

	template<int64_t i, typename E = void>
	struct sh_coeff_helper {};

	template<int64_t i>
	struct sh_coeff_helper<i, typename std::enable_if<(i & 1) == 0>::type> {
		using type = rational64<0, 1>;
	};

	template<int64_t i>
	struct sh_coeff_helper<i, typename std::enable_if<(i & 1) == 1>::type> {
		using type = rational64<1, integers::factorial<i>::val>;
	};

	template<int64_t i>
	struct sh_coeff {
		using type = typename sh_coeff_helper<i>::type;
	};

	template<int64_t i, typename E = void>
	struct cos_coeff_helper {};

	template<int64_t i>
	struct cos_coeff_helper<i, typename std::enable_if<(i & 1) == 1>::type> {
		using type = rational64<0, 1>;
	};

	template<int64_t i>
	struct cos_coeff_helper<i, typename std::enable_if<(i & 1) == 0>::type> {
		using type = rational64<integers::alternate<i / 2>::val, integers::factorial<i>::val>;
	};

	template<int64_t i>
	struct cos_coeff {
		using type = typename cos_coeff_helper<i>::type;
	};

	template<int64_t i, typename E = void>
	struct cosh_coeff_helper {};

	template<int64_t i>
	struct cosh_coeff_helper<i, typename std::enable_if<(i & 1) == 1>::type> {
		using type = rational64<0, 1>;
	};

	template<int64_t i>
	struct cosh_coeff_helper<i, typename std::enable_if<(i & 1) == 0>::type> {
		using type = rational64<1, integers::factorial<i>::val>;
	};

	template<int64_t i>
	struct cosh_coeff {
		using type = typename cosh_coeff_helper<i>::type;
	};

	template<int64_t i>
	struct geom_coeff { using type = rational64<1, 1>; };

	template<int64_t i>
	struct exp_coeff { using type = rational64<1, integers::factorial<i>::val>; };

	template<int64_t i, typename E = void>
	struct atan_coeff_helper;

	template<int64_t i>
	struct atan_coeff_helper<i, typename std::enable_if<(i & 1) == 1>::type> {
		using type = rational64<integers::alternate<i / 2>::val, i>;
	};

	template<int64_t i>
	struct atan_coeff_helper<i, typename std::enable_if<(i & 1) == 0>::type> {
		using type = rational64<0, 1>;
	};

	template<int64_t i>
	struct atan_coeff { using type = typename atan_coeff_helper<i>::type; };

	template<int64_t i, typename E = void>
	struct asin_coeff_helper;

	template<int64_t i>
	struct asin_coeff_helper<i, typename std::enable_if<(i & 1) == 1>::type>
	{
		using type = rational64<integers::factorial<i - 1>::val, integers::pow<4, i / 2>::val * integers::pow<integers::factorial<i / 2>::val, 2>::val * i>;
	};
	template<int64_t i>
	struct asin_coeff_helper<i, typename std::enable_if<(i & 1) == 0>::type>
	{
		using type = rational64<0, 1>;
	};

	template<int64_t i>
	struct asin_coeff { using type = typename asin_coeff_helper<i>::type; };

	template<int64_t i, typename E = void>
	struct asinh_coeff_helper;

	template<int64_t i>
	struct asinh_coeff_helper<i, typename std::enable_if<(i & 1) == 1>::type>
	{
		using type = rational64 < integers::alternate<i/2>::val * integers::factorial<i - 1>::val, integers::pow<4, i / 2>::val* integers::pow<integers::factorial<i / 2>::val, 2>::val* i > ;
	};
	template<int64_t i>
	struct asinh_coeff_helper<i, typename std::enable_if<(i & 1) == 0>::type>
	{
		using type = rational64<0, 1>;
	};

	template<int64_t i>
	struct asinh_coeff { using type = typename asinh_coeff_helper<i>::type; };
	
	template<int64_t i, typename E = void>
	struct atanh_coeff_helper;

	template<int64_t i>
	struct atanh_coeff_helper<i, typename std::enable_if<(i & 1) == 1>::type>
	{
		using type = rational64<1, i>;
	};
	
	template<int64_t i>
	struct atanh_coeff_helper<i, typename std::enable_if<(i & 1) == 0>::type>
	{
		using type = rational64<0, 1>;
	};

	template<int64_t i>
	struct atanh_coeff { using type = typename asinh_coeff_helper<i>::type; };

	/// ln(1+x)
	template<typename T, int64_t deg>
	using lnp1 = taylor<T, lnp1_coeff, deg>;

	/// atan(x);
	template<typename T, int64_t deg>
	using atan = taylor<T, atan_coeff, deg>;

	/// e^x
	template<typename T, int64_t deg>
	using exp = taylor<T, exp_coeff, deg>;

	/// sin(x)
	template<typename T, int64_t deg>
	using sin = taylor<T, sin_coeff, deg>;

	/// sh(x)
	template<typename T, int64_t deg>
	using sinh = taylor<T, sh_coeff, deg>;

	/// ch(x)
	template<typename T, int64_t deg>
	using cosh = taylor<T, cosh_coeff, deg>;

	/// cos(x)
	template<typename T, int64_t deg>
	using cos = taylor<T, cos_coeff, deg>;

	/// 1 / (1-x)
	template<typename T, int64_t deg>
	using geometric_sum = taylor<T, geom_coeff, deg>;

	/// asin(x)
	template<typename T, int64_t deg>
	using asin = taylor<T, asin_coeff, deg>;
	
	/// asinh(x)
	template<typename T, int64_t deg>
	using asinh = taylor<T, asinh_coeff, deg>;

	/// atanh(x)
	template<typename T, int64_t deg>
	using atanh = taylor<T, atanh_coeff, deg>;
}
