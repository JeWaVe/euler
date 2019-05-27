# euler
compile time polynomials

## supported constructs
### integer constants
- rational number
```cpp
euler::rational<p, q> r; // p/q
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

### greatest common divisor
```cpp
euler::gcd<12, 18>::val // 6
```

### polynomials
```cpp
euler::polynomial<T, rational...>
```
T can be any type castable to float and double
supports operators `+,*,-` (TODO: division)
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

### custom taylor serie
here is an example for the function [x / (e^x - 1)](https://fr.wikipedia.org/wiki/Nombre_de_Bernoulli#D%C3%A9finition_par_une_fonction_g%C3%A9n%C3%A9ratrice):

```cpp
template<int64_t i>
struct my_coeff {
	using type = decltype(euler::integers::bernouilli<i>::val / euler::rational<euler::integers::factorial<i>::val, 1>{});
};

template<typename T, int deg>
using my_func = euler::taylor<T, my_coeff, deg>;
``` 






