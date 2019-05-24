#include "euler.h"
#include <cstdio>

int main()
{
	constexpr euler::exp<double, 15> e;
	constexpr auto result = e(3.0);
	::printf("%.17lf\n", result);
}