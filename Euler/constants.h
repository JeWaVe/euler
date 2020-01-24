#pragma once
#include <cstdint>

namespace euler
{
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
}
