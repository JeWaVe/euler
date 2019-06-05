#include <algorithm>
#include <type_traits>
#include <tuple>

namespace euler {

	namespace integers {
		template<int64_t k, int64_t n>
		struct binomial;

		template<int m>
		struct bernouilli;

		// (-1)^n
		template<int64_t n>
		struct alternate;

		template<int64_t n>
		struct factorial;
	}


	template<template<auto x> typename basetype, typename a, typename b>
	struct gcd {
		static constexpr auto val = gcd<basetype, b, basetype<a::val - (a::val / b::val) * b::val>>::val;
	};

	template<template<auto x> typename basetype, typename a, typename E = void>
	struct gcd_reduce_helper;


	template<template<auto x> typename basetype, typename a>
	struct gcd_reduce_helper<basetype, a, typename std::enable_if<(a::val > 0)>::type>
	{
		static constexpr auto val = a{};
	};

	template<template<auto x> typename basetype, typename a>
	struct gcd_reduce_helper<basetype, a, typename std::enable_if<(a::val <= 0)>::type>
	{
		static constexpr auto val = -a{};
	};


	template<template<auto x> typename basetype, typename a>
	struct gcd<basetype, a, typename a::zero_type> {
		static constexpr auto val = gcd_reduce_helper<basetype, a>::val;
	};

	/* TODO
	template<typename PA, typename PB>
	struct poly_gcd {
		static constexpr auto val = poly_gcd<PB, decltype(PA{} -(PA{} / PB{}) * PB {}) > ::val;
	};

	template<typename PA, typename E = void>
	struct poly_gcd_reduce_helper;


	template<typename PA>
	struct poly_gcd_reduce_helper < PA, typename std::enable_if < (PA{} > 0) > ::type >
	{
		static constexpr auto val = PA{};
	};

	template<typename PA>
	struct poly_gcd_reduce_helper < PA, typename std::enable_if < (PA{} <= 0) > ::type >
	{
		static constexpr auto val = -PA{};
	};

	template<typename PA>
	struct poly_gcd<PA, typename PA::zero_type> {
		static constexpr auto val = poly_gcd_reduce_helper<PA>::val;
	};
	*/

	// scalar rational
	template<template<auto x> typename basetype, typename p, typename q>
	struct rational
	{
		using zero_type = rational<basetype, basetype<p::zero_type::val>, basetype<p::one_type::val>>;
		using one_type = rational<basetype, basetype<p::one_type::val>, basetype<p::one_type::val>>;

		using P = p;
		using Q = q;

		explicit constexpr operator float() const {
			return (float)p::val / (float)q::val;
		}

		explicit constexpr operator double() const {
			return (double)p::val / (double)q::val;
		}

		constexpr auto inv() const {
			return rational<basetype, q, p> {};
		}

		template<typename p1, typename q1>
		constexpr auto operator * (const rational<basetype, p1, q1>&) const {
			constexpr auto p2 = p{} *p1{};
			constexpr auto q2 = q{} *q1{};
			constexpr auto div = gcd<basetype, decltype(p2), decltype(q2)>::val;

			return rational<basetype, decltype(p2 / div), decltype(q2 / div)> {};
		}

		template<typename p1, typename q1>
		constexpr auto operator / (const rational<basetype, p1, q1>&) const {
			constexpr auto p2 = p{} *q1{};
			constexpr auto q2 = q{} *p1{};
			constexpr auto div = gcd<basetype, decltype(p2), decltype(q2)>::val;

			return rational<basetype, decltype(p2 / div), decltype(q2 / div)> {};
		}

		template<typename p1, typename q1>
		constexpr auto operator + (const rational<basetype, p1, q1>&) const {
			constexpr auto p2 = p{} *q1{} +q{} *p1{};
			constexpr auto q2 = q{} *q1{};
			constexpr auto div = gcd<basetype, decltype(p2), decltype(q2)>::val;

			return rational<basetype, decltype(p2 / div), decltype(q2 / div)> {};
		}

		template<typename p1, typename q1>
		constexpr auto operator - (const rational<basetype, p1, q1>&) const {
			constexpr auto p2 = p{} *q1{} -q{} *p1{};
			constexpr auto q2 = q{} *q1{};
			constexpr auto div = gcd<basetype, decltype(p2), decltype(q2)>::val;

			return rational<basetype, decltype(p2 / div), decltype(q2 / div)> {};
		}
	};


	template <int64_t i, typename T, typename... Ts>
	struct type_at
	{
		static_assert(i < sizeof...(Ts) + 1, "index out of range");
		using type = typename type_at<i - 1, Ts...>::type;
	};

	template <typename T, typename... Ts> struct type_at<0, T, Ts...> {
		using type = T;
	};

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

	template<typename T, template<auto x, auto y> typename rational_type, typename ...coefficients>
	struct polynomial;

	template<typename T, template<auto x, auto y> typename rational_type, typename ...As, typename ...Bs>
	constexpr auto add(polynomial<T, rational_type, As...> a, polynomial<T, rational_type, Bs...> b);

	template<typename T, template<auto x, auto y> typename rational_type, typename ...As, typename ...Bs>
	constexpr auto sub(polynomial<T, rational_type, As...> a, polynomial<T, rational_type, Bs...> b);

	template<typename T, template<auto x, auto y> typename rational_type, typename ...As, typename ...Bs>
	constexpr auto mul(polynomial<T, rational_type, As...> a, polynomial<T, rational_type, Bs...> b);

	template<typename T, template<auto x, auto y> typename rational_type, typename ...As, typename ...Bs>
	constexpr auto div(polynomial<T, rational_type, As...> a, polynomial<T, rational_type, Bs...> b);

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

	template<template<auto x, auto y> typename rational_type, int64_t N, typename ...coefficients>
	struct degree_helper {
		static constexpr auto value() {
			if constexpr (N <= 0 || N >= sizeof...(coefficients))
			{
				return 0;
			}
			else if constexpr (std::is_same<typename type_at<N, coefficients...>::type, typename type_at<N, coefficients...>::type::zero_type>::value)
			{
				return degree_helper<rational_type, N - 1, coefficients...>::value();
			}
			else
			{
				return N;
			}
		}
	};

	template<typename T, template<auto, auto> typename rational_type, typename ...coefficients>
	struct truncate_helper
	{
		template<class Typelist, typename T0>
		struct prepend;

		template<typename... Ts, typename T0>
		struct prepend<polynomial<T, rational_type, Ts...>, T0> {
			using type = polynomial<T, rational_type, T0, Ts...>;
		};

		template<int stop, int i, typename... Ts>
		struct sublist_impl {
			using type = polynomial<T, rational_type>;
		};

		template<int stop, int i, typename T0, typename... Ts>
		struct sublist_impl<stop, i, T0, Ts...>
		{
		private:
			static constexpr auto get_sublist_type()
			{
				if constexpr (i < 0)
					return typename sublist_impl<stop, i + 1, Ts...>::type{};
				else if constexpr (i < stop)
					return typename prepend<typename sublist_impl<stop, i + 1, Ts...>::type, T0>::type{};
				else
					return polynomial<T, rational_type>{};
			}

		public:
			using type = decltype(get_sublist_type());
		};

		template<int stop, typename... Ts>
		struct sublist {
			using type = typename sublist_impl<stop, 0, Ts...>::type;
		};

		template<int stop, typename... Ts>
		using sublist_t = typename sublist<stop, Ts...>::type;

		template<int deg>
		static constexpr auto apply(const polynomial<T, rational_type, coefficients...> p) {
			return sublist_t<deg, coefficients...> {};
		}
	};

	template<typename T, template<auto, auto> typename rational_type, int N, typename ...coefficients>
	struct poly_simplify {
		static constexpr auto apply(const polynomial<T, rational_type, coefficients...> p)
		{
			static_assert(sizeof...(coefficients) > 0, "invalid argument in poly_simplify");
			static_assert(N < sizeof...(coefficients), "invalid argument in poly_simplify");
			if constexpr (N <= 0)
			{
				return polynomial<T, rational_type, typename type_at<0, coefficients...>::type::zero_type> {};
			}
			else if constexpr (std::is_same<typename type_at<N, coefficients...>::type, typename type_at<0, coefficients...>::type::zero_type>::value)
			{
				return poly_simplify<T, rational_type, N - 1, coefficients...>::apply(p);
			}
			else
			{
				return truncate_helper<T, rational_type, coefficients...>::template apply<N + 1>(p);
			}
		}
	};

	template<typename T, template<auto, auto> typename rational_type, typename ...coefficients>
	constexpr auto simplify(const polynomial<T, rational_type, coefficients...> p) {
		return poly_simplify<T, rational_type, sizeof...(coefficients) - 1, coefficients...>::apply(p);
	}

	template<typename T, template<auto x, auto y> typename rational_type, typename ...coefficients>
	struct polynomial
	{
		// this is bad, as x and y could be anything. 
		using zero_type = polynomial<T, rational_type, typename rational_type<0, 1>::zero_type>;
		using one_type = polynomial<T, rational_type, typename rational_type<0, 1>::one_type>;

		using X = polynomial<T, rational_type, typename rational_type<0, 1>::zero_type, typename rational_type<0, 1>::one_type>;

		template<int64_t p>
		static constexpr auto monomial() {
			if constexpr (p <= 0)
			{
				return one_type{};
			}
			else
			{
				return monomial<p - 1>() * X {};
			}
		}

		constexpr bool is_zero() const {
			return sizeof...(coefficients) == 1 && std::is_same<typename type_at<0, coefficients...>::type, typename type_at<0, coefficients...>::type::zero_type>::value;
		}

		constexpr int64_t degree() const {
			return degree_helper<rational_type, sizeof...(coefficients) - 1, coefficients...>::value();
		}

		template<int64_t index>
		constexpr auto coeff_at() const {
			if constexpr (index < 0 || index >= sizeof...(coefficients)) {
				return typename type_at<0, coefficients...>::type::zero_type{};
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
		constexpr auto operator+(const polynomial<T, rational_type, Bs...>& b) const
		{
			return simplify(add(*this, b));
		}

		template<typename ...Bs>
		constexpr auto operator-(const polynomial<T, rational_type, Bs...>& b) const
		{
			return simplify(sub(*this, b));
		}

		template<typename ...Bs>
		constexpr auto operator*(const polynomial<T, rational_type, Bs...> b) const {
			return mul(*this, b);
		}

		template<auto x, auto y>
		constexpr auto scale(const rational_type<x, y>) const {
			return mul(*this, polynomial<T, rational_type, rational_type<x, y>>{});
		}

		template<typename ...Bs>
		constexpr auto operator/(const polynomial<T, rational_type, Bs...> b) const {
			return div(*this, b);
		}
	};

	// divides p by X^deg
	template<int64_t deg>
	struct div_monomial_helper
	{
		template<typename T, template<auto, auto> typename rational_type, typename A0, typename ...As>
		static constexpr auto apply(const polynomial<T, rational_type, A0, As...>)
		{
			return typename polynomial<T, rational_type, A0>::X{} *div_monomial_helper<deg - 1>::template apply<T, rational_type, As...>(polynomial<T, rational_type, As...> {});
		}
	};


	template<>
	struct div_monomial_helper<0>
	{
		template<typename T, template<auto, auto> typename rational_type, typename A0, typename ...As>
		static constexpr auto apply(const polynomial<T, rational_type, A0, As...> p)
		{
			return p;
		}
	};

	template<typename T, template<auto x, auto y> typename rational_type, typename ...Ns>
	struct div_helper
	{
		template<typename ...Ds>
		struct inner {
			template<typename ...Qs, typename ...Rs>
			static constexpr auto apply(const polynomial<T, rational_type, Ns...> n, const polynomial<T, rational_type, Ds...> d, const polynomial<T, rational_type, Qs...> q, const polynomial<T, rational_type, Rs...> r)
			{
				if constexpr (sizeof...(Rs) < sizeof...(Ds)) {
					return q;
				}
				else if constexpr (sizeof...(Rs) <= 1 && std::is_same<typename type_at<0, Rs...>::type::P, typename type_at<0, Rs...>::type::P::zero_type>::value)
				{
					return q;
				}
				else
				{
					constexpr auto lr = typename type_at<r.degree(), Rs...>::type{};
					constexpr auto ld = typename type_at<d.degree(), Ds...>::type{};
					constexpr auto pt = decltype(n)::template monomial<r.degree() - d.degree()>().scale(lr / ld);
					return apply(n, d, q + pt, r - pt * d);
				}
			}
		};
	};

	template<typename T, template<auto, auto> typename rational_type, typename ...Ns, typename ...Ds>
	constexpr auto div(const polynomial<T, rational_type, Ns...> n, const polynomial<T, rational_type, Ds...> d) {
		constexpr auto q = typename polynomial<T, rational_type, typename type_at<0, Ns...>::type>::zero_type{};
		return div_helper<T, rational_type, Ns...>::template inner<Ds...>::apply(n, d, q, n);
	}


	template<typename T, template<auto x, auto y> typename rational_type, int64_t k, typename ...As, typename ...Bs>
	constexpr auto add_at(const polynomial<T, rational_type, As...> a, const polynomial<T, rational_type, Bs...> b) {
		return a.template coeff_at<k>() + b.template coeff_at<k>();
	}

	template<typename T, template<auto x, auto y> typename rational_type, typename ...As, typename ...Bs, size_t ...I>
	constexpr auto add_low(polynomial<T, rational_type, As...>, polynomial<T, rational_type, Bs...>, std::index_sequence<I...>) {
		return polynomial < T, rational_type, decltype(add_at<T, rational_type, I>(polynomial<T, rational_type, As...>{}, polynomial<T, rational_type, Bs...> {}))... > {};
	}

	template<typename T, template<auto x, auto y> typename rational_type, typename ...As, typename ...Bs>
	constexpr auto add(polynomial<T, rational_type, As...> a, polynomial<T, rational_type, Bs...> b) {
		return add_low(a, b, std::make_index_sequence<std::max(sizeof...(As), sizeof...(Bs))>());
	}


	template<typename T, template<auto x, auto y> typename rational_type, int64_t k, typename ...As, typename ...Bs>
	constexpr auto sub_at(const polynomial<T, rational_type, As...> a, const polynomial<T, rational_type, Bs...> b) {
		return a.template coeff_at<k>() - b.template coeff_at<k>();
	}

	template<typename T, template<auto x, auto y> typename rational_type, typename ...As, typename ...Bs, size_t ...I>
	constexpr auto sub_low(polynomial<T, rational_type, As...>, polynomial<T, rational_type, Bs...>, std::index_sequence<I...>) {
		return polynomial < T, rational_type, decltype(sub_at<T, rational_type, I>(polynomial<T, rational_type, As...>{}, polynomial<T, rational_type, Bs...> {}))... > {};
	}

	template<typename T, template<auto x, auto y> typename rational_type, typename ...As, typename ...Bs>
	constexpr auto sub(polynomial<T, rational_type, As...> a, polynomial<T, rational_type, Bs...> b) {
		return sub_low(a, b, std::make_index_sequence<std::max(sizeof...(As), sizeof...(Bs))>());
	}

	template<typename T, template<auto x, auto y> typename rational_type, int64_t k, int64_t index, int64_t stop>
	struct mul_at_loop_helper {
		template<typename... As, typename... Bs>
		static constexpr auto apply(polynomial<T, rational_type, As...> a, polynomial<T, rational_type, Bs...> b)
		{
			return a.template coeff_at<index>() * b.template coeff_at<k - index>() + mul_at_loop_helper<T, rational_type, k, index + 1, stop>::apply(a, b);
		}
	};

	template<typename T, template<auto x, auto y> typename rational_type, int64_t k, int64_t stop>
	struct mul_at_loop_helper<T, rational_type, k, stop, stop> {
		template<typename... As, typename... Bs>
		static constexpr auto apply(polynomial<T, rational_type, As...> a, polynomial<T, rational_type, Bs...> b)
		{
			return a.template coeff_at<stop>() * b.template coeff_at<0>();
		}
	};

	template<typename T, template<auto x, auto y> typename rational_type, int64_t k, typename ...As, typename ...Bs>
	constexpr auto mul_at(polynomial<T, rational_type, As...> a, polynomial<T, rational_type, Bs...> b) {
		if constexpr (k < 0 || k >= sizeof...(As) + sizeof...(Bs))
		{
			return typename type_at<0, As...>::type::zero_type{};
		}
		else
		{
			return mul_at_loop_helper<T, rational_type, k, 0, k>::apply(a, b);
		}
	}

	template<typename T, template<auto x, auto y> typename rational_type, typename ...As, typename ...Bs, std::size_t ...I>
	constexpr auto mul_low(polynomial<T, rational_type, As...>, polynomial<T, rational_type, Bs...>, std::index_sequence<I...>) {
		return polynomial < T, rational_type, decltype(mul_at<T, rational_type, I>(polynomial<T, rational_type, As...>{}, polynomial<T, rational_type, Bs...>{}))... > {};
	}

	template<typename T, template<auto x, auto y> typename rational_type, typename ...As, typename ...Bs>
	constexpr auto mul(polynomial<T, rational_type, As...> a, polynomial<T, rational_type, Bs...> b) {
		return mul_low(a, b, std::make_index_sequence<sizeof...(As) + sizeof...(Bs) - 1>());
	}

	template<int32_t p>
	struct int32 {
		using zero_type = int32<0>;
		using one_type = int32<1>;
		static constexpr int val = p;

		constexpr auto operator-()
		{
			return int32<-p> {};
		}
	};

	template<int64_t p>
	struct int64 {
		using zero_type = int64<0>;
		using one_type = int64<1>;
		static constexpr int64_t val = p;

		constexpr auto operator-()
		{
			return int64<-p> {};
		}
	};

	template<int32_t p, int32_t q>
	constexpr auto operator*(const int32<p>&, const int32<q>&)
	{
		return int32<p * q>{};
	}

	template<int32_t p, int32_t q>
	constexpr auto operator+(const int32<p>&, const int32<q>&)
	{
		return int32<p + q>{};
	}

	template<int32_t p, int32_t q>
	constexpr auto operator-(const int32<p>&, const int32<q>&)
	{
		return int32<p - q>{};
	}

	template<int32_t p, int32_t q>
	constexpr auto operator/(const int32<p>&, const int32<q>&)
	{
		return int32<p / q>{};
	}

	template<int64_t p, int64_t q>
	constexpr auto operator*(const int64<p>&, const int64<q>&)
	{
		return int64<p * q>{};
	}

	template<int64_t p, int64_t q>
	constexpr auto operator+(const int64<p>&, const int64<q>&)
	{
		return int64<p + q>{};
	}

	template<int64_t p, int64_t q>
	constexpr auto operator-(const int64<p>&, const int64<q>&)
	{
		return int64<p - q>{};
	}

	template<int64_t p, int64_t q>
	constexpr auto operator/(const int64<p>&, const int64<q>&)
	{
		return int64<p / q>{};
	}

	template<int p, int q>
	using rational32 = rational<int32, int32<p>, int32<q>>;
	template<int64_t p, int64_t q>
	using rational64 = rational<int64, int64<p>, int64<q>>;

	template<int64_t p>
	struct integer : rational64<p, 1>
	{
		static constexpr int64_t val = p;
	};


	template<template<auto index> typename, typename, class>
	struct make_taylor_impl;

	template<template<auto index> typename coeff_at, typename T, int64_t... Is>
	struct make_taylor_impl<coeff_at, T, std::integer_sequence<int64_t, Is...>> {
		using type = polynomial<T, rational64, typename coeff_at<Is>::type...>;
	};

	// generic taylor serie, depending on coefficients
	// TODO: use that for all series below
	template<typename T, template<auto index> typename coeff_at, int64_t deg>
	using taylor = typename make_taylor_impl<coeff_at, T, std::make_integer_sequence<int64_t, deg>>::type;


	template<int64_t i>
	struct lnp1_coeff { using type = rational64<integers::alternate<i + 1>::val, i>; };

	template<>
	struct lnp1_coeff<0> { using type = rational64<0, 1>; };


	template<int64_t i, typename E = void>
	struct sin_coeff_helper {};

	template<int64_t i>
	struct sin_coeff_helper<i, typename std::enable_if<(i & 1) == 0>::type> {
		using type = rational64<0, 1>;
	};

	template<int64_t i>
	struct sin_coeff_helper<i, typename std::enable_if<(i & 1) == 1>::type> {
		using type = rational64<integers::alternate<i / 2>::val, integers::factorial<i>::val>;
	};

	template<int64_t i>
	struct sin_coeff {
		using type = typename sin_coeff_helper<i>::type;
	};

	template<int64_t i, typename E = void>
	struct sh_coeff_helper {};

	template<int64_t i>
	struct sh_coeff_helper<i, typename std::enable_if<(i & 1) == 0>::type> {
		using type = rational64<0, 1>;
	};

	template<int64_t i>
	struct sh_coeff_helper<i, typename std::enable_if<(i & 1) == 1>::type> {
		using type = rational64<1, integers::factorial<i>::val>;
	};

	template<int64_t i>
	struct sh_coeff {
		using type = typename sh_coeff_helper<i>::type;
	};

	template<int64_t i, typename E = void>
	struct cos_coeff_helper {};

	template<int64_t i>
	struct cos_coeff_helper<i, typename std::enable_if<(i & 1) == 1>::type> {
		using type = rational64<0, 1>;
	};

	template<int64_t i>
	struct cos_coeff_helper<i, typename std::enable_if<(i & 1) == 0>::type> {
		using type = rational64<integers::alternate<i / 2>::val, integers::factorial<i>::val>;
	};

	template<int64_t i>
	struct cos_coeff {
		using type = typename cos_coeff_helper<i>::type;
	};

	template<int64_t i, typename E = void>
	struct cosh_coeff_helper {};

	template<int64_t i>
	struct cosh_coeff_helper<i, typename std::enable_if<(i & 1) == 1>::type> {
		using type = rational64<0, 1>;
	};

	template<int64_t i>
	struct cosh_coeff_helper<i, typename std::enable_if<(i & 1) == 0>::type> {
		using type = rational64<1, integers::factorial<i>::val>;
	};

	template<int64_t i>
	struct cosh_coeff {
		using type = typename cosh_coeff_helper<i>::type;
	};

	template<int64_t i>
	struct geom_coeff { using type = rational64<1, 1>; };

	template<int64_t i>
	struct exp_coeff { using type = rational64<1, integers::factorial<i>::val>; };

	template<int64_t i, typename E = void>
	struct atan_coeff_helper;

	template<int64_t i>
	struct atan_coeff_helper<i, typename std::enable_if<(i & 1) == 1>::type> {
		using type = rational64<integers::alternate<i / 2>::val, i>;
	};

	template<int64_t i>
	struct atan_coeff_helper<i, typename std::enable_if<(i & 1) == 0>::type> {
		using type = rational64<0, 1>;
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
		template<int64_t k, int64_t n, typename E = void>
		struct binomial_helper;

		template<int64_t k, int64_t n>
		struct binomial_helper<k, n, typename std::enable_if<(n >= 0 && k <= n / 2 && k > 0)>::type>
		{
			static constexpr auto val = binomial_helper<k - 1, n - 1>::val * rational64<n, k>{};
		};

		template<int64_t k, int64_t n>
		struct binomial_helper<k, n, typename std::enable_if<(n >= 0 && k > (n / 2) && k > 0)>::type>
		{
			static constexpr auto val = binomial_helper<n - k, n>::val;
		};

		template<int64_t n>
		struct binomial_helper<0, n, typename std::enable_if<(n >= 0)>::type>
		{
			static constexpr auto val = rational64<1, 1>{};
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
				return bernouilli_helper<k + 1, m>::template compute(accum + binomial<k, m + 1>::val * bernouilli<k>::val);
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
			static constexpr auto val = bernouilli_helper<0, m>::template compute(rational64<0, 1>{})* rational64<-1, m + 1> {};
		};

		template<>
		struct bernouilli<0> {
			static constexpr auto val = rational64<1, 1>{};
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
}
