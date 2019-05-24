#include "euler.h"
#include <cstdio>

int main() {
	constexpr euler::exp<double, 13> e;
	constexpr auto result = e(0.1);

	::printf("%.17lf\n", result);
}