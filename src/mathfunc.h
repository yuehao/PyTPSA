//
// Created by Yue Hao on 4/4/17.
//

#ifndef PYTPSA_MATHFUNC_H
#define PYTPSA_MATHFUNC_H
#include <cstdlib>
unsigned long doublefactorial(const int&);
unsigned long factorial(const int&);
unsigned long binomial(const int& n, const int& m);

template<class T>
T gcd(const T& a, const T& b) {
    return b == 0 ? a : gcd(b, a % b);
}

template<class T>
void fraction_reduction(T& a, T& b) {
    if (a == 0 || b == 0) return;
    T temp = gcd(std::abs(a), std::abs(b));
    
    a = a / temp;
    b = b / temp;
    
    return;
}
#endif //PYTPSA_MATHFUNC_H
