#include <utility>
#include <cstdint>
#include <type_traits>
#include <algorithm>

namespace euler
{
	namespace integers {
		template<int64_t k, int64_t n>
		struct binomial;

		template <typename T, int64_t ...coefficients>
		struct polynomial;

		template<int m>
		struct bernouilli;

		// (-1)^n
		template<int64_t n>
		struct alternate;

		template<int64_t n>
		struct factorial;
	}

	template<int64_t p, int64_t q, typename E = void>
	struct rational;

	template<int64_t p>
	struct integer;

	template<int64_t a, int64_t b>
	struct gcd {
		static constexpr int64_t val = gcd<b, a - (a / b) * b>::val;
	};

	template<int64_t a>
	struct gcd<a, 0> {
		static constexpr int64_t val = a > 0 ? a : -a;
	};

	template<int64_t p, int64_t q>
	struct rational<p, q, typename ::std::enable_if<q != 0>::type> {

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

		template<int64_t p1, int64_t q1>
		constexpr auto operator * (const rational<p1, q1>&) const {
			constexpr auto p2 = p * p1;
			constexpr auto q2 = q * q1;
			constexpr auto div = gcd<p2, q2>::val;

			return rational<p2 / div, q2 / div> {};
		}

		template<int64_t p1, int64_t q1>
		constexpr auto operator / (const rational<p1, q1>&) const {
			constexpr auto p2 = p * q1;
			constexpr auto q2 = q * p1;
			constexpr auto div = gcd<p2, q2>::val;

			return rational<p2 / div, q2 / div> {};
		}

		template<int64_t p1, int64_t q1>
		constexpr auto operator + (const rational<p1, q1>&) const {
			constexpr auto p2 = p * q1 + q * p1;
			constexpr auto q2 = q * q1;
			constexpr auto div = gcd<p2, q2>::val;

			return rational<p2 / div, q2 / div> {};
		}

		template<int64_t p1, int64_t q1>
		constexpr auto operator - (const rational<p1, q1>&) const {
			constexpr auto p2 = p * q1 - q * p1;
			constexpr auto q2 = q * q1;
			constexpr auto div = gcd<p2, q2>::val;

			return rational<p2 / div, q2 / div> {};
		}
	};

	template<int64_t p>
	struct integer : rational<p, 1>
	{
		static constexpr int64_t val = p;
	};


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

	template <int64_t i, typename T, typename... Ts>
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

	template<int64_t index, int64_t stop>
	struct poly_eval_helper
	{
		template<typename T, typename ...coefficients>
		static constexpr T apply(const T& accum, const T& x)
		{
			constexpr int64_t N = sizeof...(coefficients);
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

	template<int64_t stop>
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
		template<int64_t index>
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

	template<typename T, int64_t k, typename ...As, typename ...Bs>
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

	template<typename T, int64_t k, typename ...As, typename ...Bs>
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

	template<typename T, int64_t k, int64_t index, int64_t stop>
	struct mul_at_loop_helper {
		template<typename... As, typename... Bs>
		static constexpr auto apply(polynomial<T, As...> a, polynomial<T, Bs...> b)
		{
			return a.template coeff_at<index>() * b.template coeff_at<k - index>() + mul_at_loop_helper<T, k, index + 1, stop>::apply(a, b);
		}
	};

	template<typename T, int64_t k, int64_t stop>
	struct mul_at_loop_helper<T, k, stop, stop> {
		template<typename... As, typename... Bs>
		static constexpr auto apply(polynomial<T, As...> a, polynomial<T, Bs...> b)
		{
			return a.template coeff_at<stop>() * b.template coeff_at<0>();
		}
	};

	template<typename T, int64_t k, typename ...As, typename ...Bs>
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

	template<typename T, typename ...As, typename ...Bs, std::size_t ...I>
	constexpr auto mul_low(polynomial<T, As...>, polynomial<T, Bs...>, std::index_sequence<I...>) {
		return polynomial < T, decltype(mul_at2<T, I>(polynomial<T, As...>{}, polynomial<T, Bs...>{}))... > {};
	}

	template<typename T, typename ...As, typename ...Bs>
	constexpr auto mul(polynomial<T, As...> a, polynomial<T, Bs...> b) {
		return mul_low(a, b, std::make_index_sequence<sizeof...(As) + sizeof...(Bs) - 1>());
	}
	// how to implement taylor series -- some examples

	template<template<auto index> typename, typename, class>
	struct make_taylor_impl;

	template<template<auto index> typename coeff_at, typename T, int64_t... Is>
	struct make_taylor_impl<coeff_at, T, std::integer_sequence<int64_t, Is...>> {
		using type = polynomial<T, typename coeff_at<Is>::type...>;
	};

	// generic taylor serie, depending on coefficients
	// TODO: use that for all series below
	template<typename T, template<auto index> typename coeff_at, int64_t deg>
	using taylor = typename make_taylor_impl<coeff_at, T, std::make_integer_sequence<int64_t, deg>>::type;

	template<int64_t i>
	struct lnp1_coeff { using type = rational<integers::alternate<i + 1>::val, i>; };

	template<>
	struct lnp1_coeff<0> { using type = rational<0, 1>; };


	template<int64_t i, typename E = void>
	struct sin_coeff_helper {};

	template<int64_t i>
	struct sin_coeff_helper<i, typename std::enable_if<(i & 1) == 0>::type> {
		using type = rational<0, 1>;
	};

	template<int64_t i>
	struct sin_coeff_helper<i, typename std::enable_if<(i & 1) == 1>::type> {
		using type = rational<integers::alternate<i / 2>::val, integers::factorial<i>::val>;
	};

	template<int64_t i>
	struct sin_coeff {
		using type = typename sin_coeff_helper<i>::type;
	};

	template<int64_t i, typename E = void>
	struct sh_coeff_helper {};

	template<int64_t i>
	struct sh_coeff_helper<i, typename std::enable_if<(i & 1) == 0>::type> {
		using type = rational<0, 1>;
	};

	template<int64_t i>
	struct sh_coeff_helper<i, typename std::enable_if<(i & 1) == 1>::type> {
		using type = rational<1, integers::factorial<i>::val>;
	};

	template<int64_t i>
	struct sh_coeff {
		using type = typename sh_coeff_helper<i>::type;
	};

	template<int64_t i, typename E = void>
	struct cos_coeff_helper {};

	template<int64_t i>
	struct cos_coeff_helper<i, typename std::enable_if<(i & 1) == 1>::type> {
		using type = rational<0, 1>;
	};

	template<int64_t i>
	struct cos_coeff_helper<i, typename std::enable_if<(i & 1) == 0>::type> {
		using type = rational<integers::alternate<i / 2>::val, integers::factorial<i>::val>;
	};

	template<int64_t i>
	struct cos_coeff {
		using type = typename cos_coeff_helper<i>::type;
	};

	template<int64_t i, typename E = void>
	struct cosh_coeff_helper {};

	template<int64_t i>
	struct cosh_coeff_helper<i, typename std::enable_if<(i & 1) == 1>::type> {
		using type = rational<0, 1>;
	};

	template<int64_t i>
	struct cosh_coeff_helper<i, typename std::enable_if<(i & 1) == 0>::type> {
		using type = rational<1, integers::factorial<i>::val>;
	};

	template<int64_t i>
	struct cosh_coeff {
		using type = typename cosh_coeff_helper<i>::type;
	};

	template<int64_t i>
	struct geom_coeff { using type = rational<1, 1>; };

	template<int64_t i>
	struct exp_coeff { using type = rational<1, integers::factorial<i>::val>; };

	template<int64_t i, typename E = void>
	struct atan_coeff_helper;

	template<int64_t i>
	struct atan_coeff_helper<i, typename std::enable_if<(i & 1) == 1>::type> {
		using type = rational<integers::alternate<i / 2>::val, i>;
	};

	template<int64_t i>
	struct atan_coeff_helper<i, typename std::enable_if<(i & 1) == 0>::type> {
		using type = rational<0, 1>;
	};

	template<int64_t i>
	struct atan_coeff { using type = typename atan_coeff_helper<i>::type; };

	/// ln(1+x)
	template<typename T, int64_t deg>
	using lnp1 = taylor<T, lnp1_coeff, deg>;

	/// atan(x);
	template<typename T, int64_t deg>
	using atan = taylor<T, atan_coeff, deg>;

	/// e^x
	template<typename T, int64_t deg>
	using exp = taylor<T, exp_coeff, deg>;

	/// sin(x)
	template<typename T, int64_t deg>
	using sin = taylor<T, sin_coeff, deg>;

	/// sh(x)
	template<typename T, int64_t deg>
	using sinh = taylor<T, sh_coeff, deg>;

	/// ch(x)
	template<typename T, int64_t deg>
	using cosh = taylor<T, cosh_coeff, deg>;

	/// cos(x)
	template<typename T, int64_t deg>
	using cos = taylor<T, cos_coeff, deg>;

	/// 1 / (1-x)
	template<typename T, int64_t deg>
	using geometric_sum = taylor<T, geom_coeff, deg>;

	namespace integers
	{
		template <typename T, int64_t ...coefficients>
		struct polynomial : ::euler::polynomial<T, integer<coefficients>...> { };

		template<int64_t k, int64_t n, typename E = void>
		struct binomial_helper;

		template<int64_t k, int64_t n>
		struct binomial_helper<k, n, typename std::enable_if<(n >= 0 && k <= n / 2 && k > 0)>::type> 
		{
			static constexpr auto val = binomial_helper<k - 1, n - 1>::val * rational<n, k>{};
		};

		template<int64_t k, int64_t n>
		struct binomial_helper<k, n, typename std::enable_if<(n >= 0 && k > (n / 2) && k > 0)>::type>
		{
			static constexpr auto val = binomial_helper<n - k, n>::val;
		};

		template<int64_t n>
		struct binomial_helper<0, n, typename std::enable_if<(n >= 0)>::type>
		{
			static constexpr auto val = rational<1, 1>{};
		};

		template<int64_t k, int64_t n>
		struct binomial {
			static constexpr auto val = binomial_helper<k, n>::val;
		};

		template<int m>
		struct bernouilli;

		template<int k, int m>
		struct bernouilli_helper {
			template<typename Taccum>
			static constexpr auto compute(const Taccum& accum) {
					return bernouilli_helper<k+1, m>::template compute(accum + binomial<k, m+1>::val * bernouilli<k>::val);
			}
		};

		template<int stop>
		struct bernouilli_helper<stop, stop> {
			template<typename Taccum>
			static constexpr auto compute(const Taccum& accum) {
				return accum;
			}
		};

		template<int m>
		struct bernouilli {
			static constexpr auto val = bernouilli_helper<0, m>::template compute(rational<0, 1>{}) * rational<-1, m + 1> {};
		};

		template<>
		struct bernouilli<0> {
			static constexpr auto val = rational<1, 1>{};
		};

		// (-1)^n
		template<int64_t n>
		struct alternate {
			static constexpr int64_t val = n & 1 ? -1 : 1;
		};

		template<int64_t n>
		struct factorial { static constexpr int64_t val = factorial<n - 1>::val * n; };

		template<>
		struct factorial<0> { static constexpr int64_t val = 1; };
	}
} // end namespace euler
