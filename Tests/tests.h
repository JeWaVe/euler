#pragma once
#include <algorithm>
#include <cmath>
#include <iostream>

#include "euler.h"

int test_exp() {
	constexpr double a = 0.1;
	euler::exp<double, 15> e;
	if (::std::abs((double)e(a) - ::exp(a)) > std::numeric_limits<double>::epsilon()) {
		::printf("ERROR EXP\n");
		return 1;
	}

	return 0;
}

int test_sin() {
	constexpr double a = 0.1;
	euler::sin<double, 15> e;
	if (::std::abs((double)e(a) - ::sin(a)) > std::numeric_limits<double>::epsilon()) {
		::printf("ERROR SIN\n");
		return 1;
	}

	return 0;
}

int test_cos() {
	constexpr double a = 0.1;
	euler::cos<double, 15> e;
	if (::std::abs((double)e(a) - ::cos(a)) > std::numeric_limits<double>::epsilon()) {
		::printf("ERROR COS\n");
		return 1;
	}

	return 0;
}

int test_cosh() {
	constexpr double a = 0.1;
	euler::cosh<double, 15> e;
	if (::std::abs((double)e(a) - ::cosh(a)) > std::numeric_limits<double>::epsilon()) {
		::printf("ERROR COSH\n");
		return 1;
	}

	return 0;
}

int test_sinh() {
	constexpr double a = 0.1;
	euler::sinh<double, 15> e;
	if (::std::abs((double)e(a) - ::sinh(a)) > std::numeric_limits<double>::epsilon()) {
		::printf("ERROR COSH\n");
		return 1;
	}

	return 0;
}

int test_atan() {
	constexpr double a = 0.1;
	euler::atan<double, 15> e;
	if (::std::abs((double)e(a) - ::atan(a)) > std::numeric_limits<double>::epsilon()) {
		::printf("ERROR ATAN\n");
		return 1;
	}

	return 0;
}


