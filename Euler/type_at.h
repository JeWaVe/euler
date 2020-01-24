#pragma once
#include <cstdint>

namespace euler
{
	template <int64_t i, typename T, typename... Ts>
	struct type_at
	{
		static_assert(i < sizeof...(Ts) + 1, "index out of range");
		using type = typename type_at<i - 1, Ts...>::type;
	};

	template <typename T, typename... Ts> struct type_at<0, T, Ts...> {
		using type = T;
	};
}
