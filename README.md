# deprecated
This repo is NOT maintained anymore - all code has moved to [aerobus](https://github.com/aerobus-open-source/aerobus)

# euler
compile time polynomials

## requirements
- this is a single-header file (euler.h)
- any c++-17 compatible compiler

## supported constructs
### integer constants
- rational number
```cpp
euler::rational64<p, q> r; // p/q
```
supports following operators : `*,+,/,-`
- binomial coefficients
```cpp
euler::integers::binomial<k, n>::val // C(k,n) (rational)
```
- bernouilli numbers
```cpp
euler::integers::bernouilli<m>::val // B_m (rational)
```
- alternate
```cpp
euler::integers::alternate<k>::val // (-1)^k (int64_t)
```
- factorial
```cpp
euler::integers::factorial<n>::val // n! (int64_t)
```

- natural powers
```cpp
euler::integers::pow<n, p>::val // n^p

### greatest common divisor
```cpp
euler::gcd<euler::int32<12>, euler::int32<18>>::val // 6
```

### polynomials
```cpp
euler::polynomial64<T, rational...>
```
T can be any type castable to float and double
supports operators `+,*,-,/`
supports evaluation on compile time constant 

### taylor series
Some important functions are expressed as polynomial (taylor series in zero), such as
- exp
- lnp1 (ln(1+x))
- sin
- cos
- atan
- sinh
- cosh
- geometric_sum (1/(1-x))
- asin
- asinh
- atanh

### custom taylor serie
here is an example for the function [x / (e^x - 1)](https://fr.wikipedia.org/wiki/Nombre_de_Bernoulli#D%C3%A9finition_par_une_fonction_g%C3%A9n%C3%A9ratrice):

```cpp
template<int64_t i>
struct my_coeff {
	using type = decltype(integers::bernouilli<i>::val / rational64<integers::factorial<i>::val, 1>{});
};

template<int deg>
using my_func = taylor<double, my_coeff, deg>;
``` 
Example in VS 2019:
![toto](./images/taylor.png)





