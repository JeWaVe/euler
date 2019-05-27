#include "euler.h"
#include <cstdio>

template<int64_t index>
struct my_coeff {
	using type = euler::rational<1, 1>;
};

// 1/(1-x)
template<typename T, int64_t d>
using serie = euler::taylor<T, my_coeff, d>;

int main()
{
	constexpr auto a = euler::integers::binomial<12, 128>::val;
	constexpr serie<double, 15> e;
	constexpr auto result = e(0.1);
	::printf("%.17lf\n", result);
}
