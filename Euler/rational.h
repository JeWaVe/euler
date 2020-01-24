#pragma once
#include "defines.h"
#include "forward_declarations.h"

namespace euler {
	// scalar rational
	template<template<auto x> typename basetype, typename p, typename q>
	struct rational
	{
		using zero_type = rational<basetype, basetype<p::zero_type::val>, basetype<p::one_type::val>>;
		using one_type = rational<basetype, basetype<p::one_type::val>, basetype<p::one_type::val>>;

		using P = p;
		using Q = q;

		explicit constexpr operator float() const {
			return static_cast<float>(p::val) / static_cast<float>(q::val);
		}

		explicit constexpr operator double() const {
			return static_cast<double>(p::val) / static_cast<double>(q::val);
		}

		CONST inv() const {
			return rational<basetype, q, p> {};
		}

		template<typename p1, typename q1>
		CONST operator * (const rational<basetype, p1, q1>&) const {
			CONST p2 = p{} *p1{};
			CONST q2 = q{} *q1{};
			CONST div = gcd<basetype, decltype(p2), decltype(q2)>::val;

			return rational<basetype, decltype(p2 / div), decltype(q2 / div)> {};
		}

		template<typename p1, typename q1>
		CONST operator / (const rational<basetype, p1, q1>&) const {
			CONST p2 = p{} *q1{};
			CONST q2 = q{} *p1{};
			CONST div = gcd<basetype, decltype(p2), decltype(q2)>::val;

			return rational<basetype, decltype(p2 / div), decltype(q2 / div)> {};
		}

		template<typename p1, typename q1>
		CONST operator + (const rational<basetype, p1, q1>&) const {
			CONST p2 = p{} *q1{} +q{} *p1{};
			CONST q2 = q{} *q1{};
			CONST div = gcd<basetype, decltype(p2), decltype(q2)>::val;

			return rational<basetype, decltype(p2 / div), decltype(q2 / div)> {};
		}

		template<typename p1, typename q1>
		CONST operator - (const rational<basetype, p1, q1>&) const {
			CONST p2 = p{} *q1{} -q{} *p1{};
			CONST q2 = q{} *q1{};
			CONST div = gcd<basetype, decltype(p2), decltype(q2)>::val;

			return rational<basetype, decltype(p2 / div), decltype(q2 / div)> {};
		}
	};
}
