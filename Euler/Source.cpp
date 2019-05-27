#include "euler.h"

// how to implement a taylor serie for x / (e^x - 1) (https://fr.wikipedia.org/wiki/Nombre_de_Bernoulli#D%C3%A9finition_par_une_fonction_g%C3%A9n%C3%A9ratrice)

template<int64_t i>
struct my_coeff {
	using type = decltype(euler::integers::bernouilli<i>::val / euler::rational<euler::integers::factorial<i>::val, 1>{});
};

template<typename T, int deg>
using my_func = euler::taylor<T, my_coeff, deg>;

int main() {
	my_func<double, 15> f;
	constexpr auto a = f(0.1);
}