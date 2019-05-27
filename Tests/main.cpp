#include <cstdlib>
#include <cstdio>

#include "tests.h"


int main(int argc, char* argv[]) {
	if (test_exp() != 0)
		return 1;
	if (test_sin() != 0)
		return 1;
	if (test_cos() != 0)
		return 1;
	if (test_sinh() != 0)
		return 1;
	if (test_cosh() != 0)
		return 1;
	if (test_atan() != 0)
		return 1;
	::printf("OK\n");
	return 0;
}