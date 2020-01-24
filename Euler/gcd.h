#pragma once
#include <type_traits>
#include "defines.h"

namespace euler
{
	template<template<auto x> typename basetype, typename a, typename b>
	struct gcd {
		static CONST val = gcd<basetype, b, basetype<a::val - (a::val / b::val) * b::val>>::val;
	};

	template<template<auto x> typename basetype, typename a, typename E = void>
	struct gcd_reduce_helper;


	template<template<auto x> typename basetype, typename a>
	struct gcd_reduce_helper<basetype, a, typename std::enable_if<(a::val > 0)>::type>
	{
		static CONST val = a{};
	};

	template<template<auto x> typename basetype, typename a>
	struct gcd_reduce_helper<basetype, a, typename std::enable_if<(a::val <= 0)>::type>
	{
		static CONST val = -a{};
	};


	template<template<auto x> typename basetype, typename a>
	struct gcd<basetype, a, typename a::zero_type> {
		static CONST val = gcd_reduce_helper<basetype, a>::val;
	};

	/* TODO
	template<typename PA, typename PB>
	struct poly_gcd {
		static CONST val = poly_gcd<PB, decltype(PA{} -(PA{} / PB{}) * PB {}) > ::val;
	};

	template<typename PA, typename E = void>
	struct poly_gcd_reduce_helper;


	template<typename PA>
	struct poly_gcd_reduce_helper < PA, typename std::enable_if < (PA{} > 0) > ::type >
	{
		static CONST val = PA{};
	};

	template<typename PA>
	struct poly_gcd_reduce_helper < PA, typename std::enable_if < (PA{} <= 0) > ::type >
	{
		static CONST val = -PA{};
	};

	template<typename PA>
	struct poly_gcd<PA, typename PA::zero_type> {
		static CONST val = poly_gcd_reduce_helper<PA>::val;
	};
	*/
}
