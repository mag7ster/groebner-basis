## Grobner Basis
Modern C++ library for working with polynomials of several variables and searching Groebner Basises.

What is implemented?
  1. Monomials
  2. Monomial orders -- Lex, RevLex, GrLex, GrevLex
  3. Terms (Monomial with coefficient)
  4. Polynomials
  5. Set of Polynomials
  6. Buchberger Algorithm 
 
# Build

Install vcpkg before building
 ```bash
mkdir build
cd build
cmake --preset=default ..
```

# Compile and run benchmark
```bash
make bench
./bench
```


# Compile and run tests
```bash
python3 ../generate_tests.py
make tests
./tests
```
You may need to install SymPy for python to generate tests
