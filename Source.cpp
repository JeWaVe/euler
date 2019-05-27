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
	constexpr euler::atan<double, 15> e;
	constexpr auto result = e(0.1);

	constexpr auto zou = euler::integers::bernouilli<12>::val;
	::printf("%.17lf\n", result);
}
