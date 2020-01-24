#pragma once
#include <cstdint>
#include "defines.h"
#include "forward_declarations.h"

namespace euler
{
	template<int32_t p>
	struct int32 {
		using zero_type = int32<0>;
		using one_type = int32<1>;
		static constexpr int val = p;

		CONST operator-() const
		{
			return int32<-p> {};
		}
	};

	template<int64_t p>
	struct int64 {
		using zero_type = int64<0>;
		using one_type = int64<1>;
		static constexpr int64_t val = p;

		CONST operator-() const
		{
			return int64<-p> {};
		}
	};

	template<int32_t p, int32_t q>
	CONST operator*(const int32<p>&, const int32<q>&)
	{
		return int32<p * q>{};
	}

	template<int32_t p, int32_t q>
	CONST operator+(const int32<p>&, const int32<q>&)
	{
		return int32<p + q>{};
	}

	template<int32_t p, int32_t q>
	CONST operator-(const int32<p>&, const int32<q>&)
	{
		return int32<p - q>{};
	}

	template<int32_t p, int32_t q>
	CONST operator/(const int32<p>&, const int32<q>&)
	{
		return int32<p / q>{};
	}

	template<int64_t p, int64_t q>
	CONST operator*(const int64<p>&, const int64<q>&)
	{
		return int64<p * q>{};
	}

	template<int64_t p, int64_t q>
	CONST operator+(const int64<p>&, const int64<q>&)
	{
		return int64<p + q>{};
	}

	template<int64_t p, int64_t q>
	CONST operator-(const int64<p>&, const int64<q>&)
	{
		return int64<p - q>{};
	}

	template<int64_t p, int64_t q>
	CONST operator/(const int64<p>&, const int64<q>&)
	{
		return int64<p / q>{};
	}

	template<int p, int q>
	using rational32 = rational<int32, int32<p>, int32<q>>;
	template<int64_t p, int64_t q>
	using rational64 = rational<int64, int64<p>, int64<q>>;

	template<typename T, typename ...coefficients>
	using polynomial32 = polynomial<T, rational32, coefficients...>;

	template<typename T, typename ...coefficients>
	using polynomial64 = polynomial<T, rational64, coefficients...>;

	template<int64_t p>
	struct integer : rational64<p, 1>
	{
		static constexpr int64_t val = p;
	};

	namespace integers
	{
		template<int64_t k, int64_t n, typename E = void>
		struct binomial_helper;

		template<int64_t k, int64_t n>
		struct binomial_helper<k, n, typename std::enable_if<(n >= 0 && k <= n / 2 && k > 0)>::type>
		{
			static CONST val = binomial_helper<k - 1, n - 1>::val * rational64<n, k>{};
		};

		template<int64_t k, int64_t n>
		struct binomial_helper<k, n, typename std::enable_if<(n >= 0 && k > (n / 2) && k > 0)>::type>
		{
			static CONST val = binomial_helper<n - k, n>::val;
		};

		template<int64_t n>
		struct binomial_helper<0, n, typename std::enable_if<(n >= 0)>::type>
		{
			static CONST val = rational64<1, 1>{};
		};

		template<int64_t k, int64_t n>
		struct binomial {
			static CONST val = binomial_helper<k, n>::val;
		};

		template<int m>
		struct bernouilli;

		template<int k, int m>
		struct bernouilli_helper {
			template<typename Taccum>
			static CONST compute(const Taccum& accum) {
				return bernouilli_helper<k + 1, m>::template compute(accum + binomial<k, m + 1>::val * bernouilli<k>::val);
			}
		};

		template<int stop>
		struct bernouilli_helper<stop, stop> {
			template<typename Taccum>
			static CONST compute(const Taccum& accum) {
				return accum;
			}
		};

		template<int m>
		struct bernouilli {
			static CONST val = bernouilli_helper<0, m>::template compute(rational64<0, 1>{})* rational64<-1, m + 1> {};
		};

		template<>
		struct bernouilli<0> {
			static CONST val = rational64<1, 1>{};
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

		template<int64_t p, int64_t n>
		struct pow { static constexpr int64_t val = p * pow<p, n - 1>::val; };

		template<int64_t p>
		struct pow<p, 0> { static constexpr int64_t val = 1; };
	}
}
