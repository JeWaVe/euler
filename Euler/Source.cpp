#include "euler.h"

using namespace euler;
// how to implement a taylor serie for x / (e^x - 1) (https://fr.wikipedia.org/wiki/Nombre_de_Bernoulli#D%C3%A9finition_par_une_fonction_g%C3%A9n%C3%A9ratrice)

template<int64_t i>
struct my_coeff {
	using type = decltype(integers::bernouilli<i>::val / rational64<integers::factorial<i>::val, 1>{});
};

template<int deg>
using my_func = taylor<double, my_coeff, deg>;

int main() {
	constexpr auto a = my_func<8>{}; // a / B leads to integers overflow when using 15 instead of 8
	constexpr auto b = polynomial64<double, rational64<-3, 1>, rational64<1, 1>>{};
	constexpr auto q = a / b; // WARNING : division is made on the polynoms here, so q(x) != a(x) / b(x)
	::printf("WARN : %lf == %lf\n", q(0.1), a(0.1) / b(0.1));

	constexpr auto r = a - b * q;
	constexpr auto aa = b * q + r; // aa == a
	::printf(" BUT : %lf == %lf\n", (aa)(0.1), a(0.1));
}