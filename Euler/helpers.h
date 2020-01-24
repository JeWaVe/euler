#pragma once
namespace euler
{
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
				constexpr T coeff = static_cast<T>(typename type_at<index, coefficients...>::type{});
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
		static CONST value() {
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
			static CONST get_sublist_type()
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
		static CONST apply(const polynomial<T, rational_type, coefficients...> p) {
			return sublist_t<deg, coefficients...> {};
		}
	};

	template<typename T, template<auto, auto> typename rational_type, int N, typename ...coefficients>
	struct poly_simplify {
		static CONST apply(const polynomial<T, rational_type, coefficients...> p)
		{
			static_assert(sizeof...(coefficients) > 0, "invalid argument in poly_simplify");
			static_assert(N < sizeof...(coefficients), "invalid argument in poly_simplify");
			if constexpr (N < 0)
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
	CONST simplify(const polynomial<T, rational_type, coefficients...> p) {
		return poly_simplify<T, rational_type, sizeof...(coefficients) - 1, coefficients...>::apply(p);
	}

	template<typename T, template<auto x, auto y> typename rational_type, typename ...Ns>
	struct div_helper
	{
		template<typename ...Ds>
		struct inner {
			template<typename ...Qs, typename ...Rs>
			static CONST apply(const polynomial<T, rational_type, Ns...> n, const polynomial<T, rational_type, Ds...> d, const polynomial<T, rational_type, Qs...> q, const polynomial<T, rational_type, Rs...> r)
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
					CONST lr = typename type_at<r.degree(), Rs...>::type{};
					CONST ld = typename type_at<d.degree(), Ds...>::type{};
					CONST pt = decltype(n)::template monomial<r.degree() - d.degree()>().scale(lr / ld);
					return apply(n, d, q + pt, r - pt * d);
				}
			}
		};
	};

	template<typename T, template<auto x, auto y> typename rational_type, int64_t k, typename ...As, typename ...Bs>
	CONST add_at(const polynomial<T, rational_type, As...> a, const polynomial<T, rational_type, Bs...> b) {
		return a.template coeff_at<k>() + b.template coeff_at<k>();
	}

	template<typename T, template<auto x, auto y> typename rational_type, typename ...As, typename ...Bs, size_t ...I>
	CONST add_low(polynomial<T, rational_type, As...>, polynomial<T, rational_type, Bs...>, std::index_sequence<I...>) {
		return polynomial < T, rational_type, decltype(add_at<T, rational_type, I>(polynomial<T, rational_type, As...>{}, polynomial<T, rational_type, Bs...> {}))... > {};
	}

	template<typename T, template<auto x, auto y> typename rational_type, int64_t k, typename ...As, typename ...Bs>
	CONST sub_at(const polynomial<T, rational_type, As...> a, const polynomial<T, rational_type, Bs...> b) {
		return a.template coeff_at<k>() - b.template coeff_at<k>();
	}

	template<typename T, template<auto x, auto y> typename rational_type, typename ...As, typename ...Bs, size_t ...I>
	CONST sub_low(polynomial<T, rational_type, As...>, polynomial<T, rational_type, Bs...>, std::index_sequence<I...>) {
		return polynomial < T, rational_type, decltype(sub_at<T, rational_type, I>(polynomial<T, rational_type, As...>{}, polynomial<T, rational_type, Bs...> {}))... > {};
	}

	template<typename T, template<auto x, auto y> typename rational_type, int64_t k, int64_t index, int64_t stop>
	struct mul_at_loop_helper {
		template<typename... As, typename... Bs>
		static CONST apply(polynomial<T, rational_type, As...> a, polynomial<T, rational_type, Bs...> b)
		{
			return a.template coeff_at<index>() * b.template coeff_at<k - index>() + mul_at_loop_helper<T, rational_type, k, index + 1, stop>::apply(a, b);
		}
	};

	template<typename T, template<auto x, auto y> typename rational_type, int64_t k, int64_t stop>
	struct mul_at_loop_helper<T, rational_type, k, stop, stop> {
		template<typename... As, typename... Bs>
		static CONST apply(polynomial<T, rational_type, As...> a, polynomial<T, rational_type, Bs...> b)
		{
			return a.template coeff_at<stop>() * b.template coeff_at<0>();
		}
	};

	template<typename T, template<auto x, auto y> typename rational_type, int64_t k, typename ...As, typename ...Bs>
	CONST mul_at(polynomial<T, rational_type, As...> a, polynomial<T, rational_type, Bs...> b) {
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
	CONST mul_low(polynomial<T, rational_type, As...>, polynomial<T, rational_type, Bs...>, std::index_sequence<I...>) {
		return polynomial < T, rational_type, decltype(mul_at<T, rational_type, I>(polynomial<T, rational_type, As...>{}, polynomial<T, rational_type, Bs...>{}))... > {};
	}

}