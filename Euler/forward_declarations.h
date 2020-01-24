#pragma once
#include <cstdint>
#include "defines.h"

namespace euler
{
	namespace integers {
		template<int64_t k, int64_t n>
		struct binomial;

		template<int m>
		struct bernouilli;

		// (-1)^n
		template<int64_t n>
		struct alternate;

		template<int64_t n>
		struct factorial;
	}

	template<template<auto x> typename basetype, typename p, typename q>
	struct rational;

	template<template<auto x> typename basetype, typename a, typename b>
	struct gcd;

	template<typename T, template<auto x, auto y> typename rational_type, typename ...coefficients>
	struct polynomial;

	template<typename T, template<auto x, auto y> typename rational_type, typename ...As, typename ...Bs>
	CONST add(polynomial<T, rational_type, As...> a, polynomial<T, rational_type, Bs...> b);

	template<typename T, template<auto x, auto y> typename rational_type, typename ...As, typename ...Bs>
	CONST sub(polynomial<T, rational_type, As...> a, polynomial<T, rational_type, Bs...> b);

	template<typename T, template<auto x, auto y> typename rational_type, typename ...As, typename ...Bs>
	CONST mul(polynomial<T, rational_type, As...> a, polynomial<T, rational_type, Bs...> b);

	template<typename T, template<auto x, auto y> typename rational_type, typename ...As, typename ...Bs>
	CONST div(polynomial<T, rational_type, As...> a, polynomial<T, rational_type, Bs...> b);

}
