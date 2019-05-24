#pragma once

#include <utility>
#include <cstdint>
#include <type_traits>

namespace euler
{
	template<int p, int q, typename E = void>
	struct rational;

	template<int p>
	struct integer;

	template<int a, int b>
	struct gcd {
		static constexpr int val = gcd<b, a - (a / b) * b>::val;
	};

	template<int a>
	struct gcd<a, 0> {
		static constexpr int val = a > 0 ? a : -a;
	};

	template<int p, int q>
	struct rational<p, q, typename std::enable_if<q != 0>::type> {

		explicit constexpr operator float() const {
			return (float)p / (float)q;
		}

		explicit constexpr operator double() const {
			return (double)p / (double)q;
		}

		constexpr float operator * (const float& x) const {
			return ((float)p * x) / (float)q;
		}

		constexpr double operator * (const double& x) const {
			return ((double)p * x) / (double)q;
		}

		template<int p1, int q1>
		constexpr auto operator * (const rational<p1, q1>&) const {
			constexpr auto p2 = p * p1;
			constexpr auto q2 = q * q1;
			constexpr auto div = gcd<p2, q2>::val;

			return rational<p2 / div, q2 / div> {};
		}

		template<int p1, int q1>
		constexpr auto operator / (const rational<p1, q1>&) const {
			constexpr auto p2 = p * q1;
			constexpr auto q2 = q * p1;
			constexpr auto div = gcd<p2, q2>::val;

			return rational<p2 / div, q2 / div> {};
		}

		template<int p1, int q1>
		constexpr auto operator + (const rational<p1, q1>&) const {
			constexpr auto p2 = p * q1 + q * p1;
			constexpr auto q2 = q * q1;
			constexpr auto div = gcd<p2, q2>::val;

			return rational<p2 / div, q2 / div> {};
		}

		template<int p1, int q1>
		constexpr auto operator - (const rational<p1, q1>&) const {
			constexpr auto p2 = p * q1 - q * p1;
			constexpr auto q2 = q * q1;
			constexpr auto div = gcd<p2, q2>::val;

			return rational<p2 / div, q2 / div> {};
		}
	};

	template<int p>
	struct integer : rational<p, 1> {};


	template<typename T, typename ...coefficients>
	struct polynomial;

	template<typename T>
	struct constants;

	template<> struct constants<float> { static constexpr float zero = 0.0F; };
	template<> struct constants<double> { static constexpr double zero = 0.0; };
	template<> struct constants<int32_t> { static constexpr int32_t zero = 0; };
	template<> struct constants<int64_t> { static constexpr int64_t zero = 0LL; };
	template<> struct constants<short> { static constexpr short zero = 0; };
	template<> struct constants<unsigned short> { static constexpr unsigned short zero = 0; };
	template<> struct constants<char> { static constexpr char zero = 0; };
	template<> struct constants<unsigned char> { static constexpr unsigned char zero = 0; };

	template <int i, typename T, typename... Ts>
	struct type_at
	{
		static_assert(i < sizeof...(Ts) + 1, "index out of range");
		typedef typename type_at<i - 1, Ts...>::type type;
	};

	template <typename T, typename... Ts> struct type_at<0, T, Ts...> {
		typedef T type;
	};

	template<typename T, typename ...As, typename ...Bs>
	constexpr auto add(polynomial<T, As...> a, polynomial<T, Bs...> b);

	template<typename T, typename ...As, typename ...Bs>
	constexpr auto sub(polynomial<T, As...> a, polynomial<T, Bs...> b);

	template<typename T, typename ...As, typename ...Bs>
	constexpr auto mul(polynomial<T, As...> a, polynomial<T, Bs...> b);

	template<int index, int stop>
	struct poly_eval_helper
	{
		template<typename T, typename ...coefficients>
		static constexpr T apply(const T& accum, const T& x)
		{
			constexpr int N = sizeof...(coefficients);
			if constexpr (N == 0)
			{
				return constants<T>::zero;
			}
			else
			{
				constexpr T coeff = (T) typename type_at<index, coefficients...>::type{};
				return poly_eval_helper<index - 1, stop>::template apply<T, coefficients...>(x * accum + coeff, x);
			}
		}
	};

	template<int stop>
	struct poly_eval_helper<-1, stop> {
		template<typename T, typename ...coefficients>
		static constexpr T apply(const T& accum, const T& x)
		{
			return accum;
		}
	};

	template<typename T, typename ...coefficients>
	struct polynomial
	{
		template<int index>
		constexpr auto coeff_at() const {
			if constexpr (index < 0 || index >= sizeof...(coefficients)) {
				return integer<0> {};
			}
			else {
				return typename type_at<index, coefficients...>::type{};
			}
		}

		constexpr T operator()(const T& x) const {
			if constexpr (sizeof...(coefficients) == 0)
			{
				return constants<T>::zero;
			}
			else
			{
				return poly_eval_helper<sizeof...(coefficients) - 1, sizeof...(coefficients)>::template apply<T, coefficients...>(constants<T>::zero, x);
			}
		}

		template<typename ...Bs>
		constexpr auto operator+(const polynomial<T, Bs...> b) const {
			return add(*this, b);
		}

		template<typename ...Bs>
		constexpr auto operator-(const polynomial<T, Bs...> b) const {
			return sub(*this, b);
		}

		template<typename ...Bs>
		constexpr auto operator*(const polynomial<T, Bs...> b) const {
			return mul(*this, b);
		}
	};

	template<typename T, int k, typename ...As, typename ...Bs>
	constexpr auto add_at(const polynomial<T, As...> a, const polynomial<T, Bs...> b) {
		return a.template coeff_at<k>() + b.template coeff_at<k>();
	}

	template<typename T, typename ...As, typename ...Bs, size_t ...I>
	constexpr auto add_low(polynomial<T, As...>, polynomial<T, Bs...>, std::index_sequence<I...>) {
		return polynomial < T, decltype(add_at<T, I>(polynomial<T, As...>{}, polynomial<T, Bs...> {}))... > {};
	}

	template<typename T, typename ...As, typename ...Bs>
	constexpr auto add(polynomial<T, As...> a, polynomial<T, Bs...> b) {
		return add_low(a, b, std::make_index_sequence<std::max(sizeof...(As), sizeof...(Bs))>());
	}

	template<typename T, int k, typename ...As, typename ...Bs>
	constexpr auto sub_at(const polynomial<T, As...> a, const polynomial<T, Bs...> b) {
		return a.template coeff_at<k>() - b.template coeff_at<k>();
	}

	template<typename T, typename ...As, typename ...Bs, size_t ...I>
	constexpr auto sub_low(polynomial<T, As...>, polynomial<T, Bs...>, std::index_sequence<I...>) {
		return polynomial < T, decltype(sub_at<T, I>(polynomial<T, As...>{}, polynomial<T, Bs...> {}))... > {};
	}

	template<typename T, typename ...As, typename ...Bs>
	constexpr auto sub(polynomial<T, As...> a, polynomial<T, Bs...> b) {
		return sub_low(a, b, std::make_index_sequence<sizeof...(As) + sizeof...(Bs) - 1>());
	}

	template<typename T, int k, int index, int stop>
	struct mul_at_loop_helper {
		template<typename... As, typename... Bs>
		static constexpr auto apply(polynomial<T, As...> a, polynomial<T, Bs...> b)
		{
			return a.template coeff_at<index>() * b.template coeff_at<k - index>() + mul_at_loop_helper<T, k, index + 1, stop>::apply(a, b);
		}
	};

	template<typename T, int k, int stop>
	struct mul_at_loop_helper<T, k, stop, stop> {
		template<typename... As, typename... Bs>
		static constexpr auto apply(polynomial<T, As...> a, polynomial<T, Bs...> b)
		{
			return a.template coeff_at<stop>() * b.template coeff_at<0>();
		}
	};

	template<typename T, int k, typename ...As, typename ...Bs>
	constexpr auto mul_at2(polynomial<T, As...> a, polynomial<T, Bs...> b) {
		if constexpr (k < 0 || k >= sizeof...(As) + sizeof...(Bs))
		{
			return integer<0>{};
		}
		else
		{
			return mul_at_loop_helper<T, k, 0, k>::apply(a, b);
		}
	}

	template<typename T, typename ...As, typename ...Bs, size_t ...I>
	constexpr auto mul_low(polynomial<T, As...>, polynomial<T, Bs...>, std::index_sequence<I...>) {
		return polynomial < T, decltype(mul_at2<T, I>(polynomial<T, As...>{}, polynomial<T, Bs...>{}))... > {};
	}

	template<typename T, typename ...As, typename ...Bs>
	constexpr auto mul(polynomial<T, As...> a, polynomial<T, Bs...> b) {
		return mul_low(a, b, std::make_index_sequence<sizeof...(As) + sizeof...(Bs) - 1>());
	}

	// (-1)^n
	template<int n>
	struct alternate {
		static constexpr int64_t val = n & 1 ? -1 : 1;
	};

	template<int64_t n>
	struct factorial { static constexpr int64_t val = factorial<n - 1>::val * n; };

	template<>
	struct factorial<0> { static constexpr int val = 1; };

	// how to implement taylor series -- some examples
	template<class, class>
	struct make_exp_impl;

	template<class, class>
	struct make_lnp1_impl;

	template<class, class>
	struct make_sin_impl;

	template<class, class>
	struct make_sh_impl;

	template<class, class>
	struct make_cos_impl;

	template<class, class>
	struct make_cosh_impl;

	template<class, class>
	struct make_geom_impl;

	template<typename T, int... Is>
	struct make_exp_impl<T, std::integer_sequence<int, Is...>> {
		using type = polynomial<T, rational<1, factorial<Is>::val>...>;
	};

	template<typename T, int... Is>
	struct make_lnp1_impl<T, std::integer_sequence<int, Is...>> {
		using type = polynomial<T, rational<0, 1>, rational<alternate<Is>::val, Is + 1>...>;
	};

	template<int i, typename E = void>
	struct sin_coeff {};

	template<int i>
	struct sin_coeff<i, typename std::enable_if<(i & 1) == 0>::type> {
		using type = rational<0, 1>;
	};

	template<int i>
	struct sin_coeff<i, typename std::enable_if<(i & 1) == 1>::type> {
		using type = rational<alternate<i / 2>::val, factorial<i>::val>;
	};

	template<typename T, int... Is>
	struct make_sin_impl<T, std::integer_sequence<int, Is...>> {
		using type = polynomial<T, typename sin_coeff<Is>::type...>;
	};

	template<int i, typename E = void>
	struct sh_coeff {};

	template<int i>
	struct sh_coeff<i, typename std::enable_if<(i & 1) == 0>::type> {
		using type = rational<0, 1>;
	};

	template<int i>
	struct sh_coeff<i, typename std::enable_if<(i & 1) == 1>::type> {
		using type = rational<1, factorial<i>::val>;
	};

	template<typename T, int... Is>
	struct make_sh_impl<T, std::integer_sequence<int, Is...>> {
		using type = polynomial<T, typename sh_coeff<Is>::type...>;
	};

	template<int i, typename E = void>
	struct cos_coeff {};

	template<int i>
	struct cos_coeff<i, typename std::enable_if<(i & 1) == 1>::type> {
		using type = rational<0, 1>;
	};

	template<int i>
	struct cos_coeff<i, typename std::enable_if<(i & 1) == 0>::type> {
		using type = rational<alternate<i / 2>::val, factorial<i>::val>;
	};

	template<typename T, int... Is>
	struct make_cos_impl<T, std::integer_sequence<int, Is...>> {
		using type = polynomial<T, typename cos_coeff<Is>::type...>;
	};


	template<int i, typename E = void>
	struct cosh_coeff {};

	template<int i>
	struct cosh_coeff<i, typename std::enable_if<(i & 1) == 1>::type> {
		using type = rational<0, 1>;
	};

	template<int i>
	struct cosh_coeff<i, typename std::enable_if<(i & 1) == 0>::type> {
		using type = rational<1, factorial<i>::val>;
	};

	template<typename T, int... Is>
	struct make_cosh_impl<T, std::integer_sequence<int, Is...>> {
		using type = polynomial<T, typename cosh_coeff<Is>::type...>;
	};

	template<int i, typename E = void>
	struct geom_coeff { using type = rational<1, 1>; };

	template<typename T, int... Is>
	struct make_geom_impl<T, std::integer_sequence<int, Is...>> {
		using type = polynomial<T, typename geom_coeff<Is>::type...>;
	};

	/// ln(1+x)
	template<typename T, int deg>
	using lnp1 = typename make_lnp1_impl<T, std::make_integer_sequence<int, deg>>::type;

	/// e^x
	template<typename T, int deg>
	using exp = typename make_exp_impl<T, std::make_integer_sequence<int, deg>>::type;

	/// sin(x)
	template<typename T, int deg>
	using sin = typename make_sin_impl<T, std::make_integer_sequence<int, deg>>::type;

	/// sh(x)
	template<typename T, int deg>
	using sinh = typename make_sh_impl<T, std::make_integer_sequence<int, deg>>::type;

	/// ch(x)
	template<typename T, int deg>
	using cosh = typename make_cosh_impl<T, std::make_integer_sequence<int, deg>>::type;

	/// cos(x)
	template<typename T, int deg>
	using cos = typename make_cos_impl<T, std::make_integer_sequence<int, deg>>::type;

	/// 1 / (1-x)
	template<typename T, int deg>
	using geometric_sum = typename make_geom_impl<T, std::make_integer_sequence<int, deg>>::type;

	namespace integers
	{
		template <typename T, int ...coefficients>
		struct polynomial : ::euler::polynomial<T, integer<coefficients>...> { };
	}
} // end namespace hybridizer