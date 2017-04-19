//
// Created by Yue Hao on 4/5/17.
//

#ifndef PYTPSA_RATIONAL_H
#define PYTPSA_RATIONAL_H
#include "mathfunc.h"
#include <ostream>

template<class T>
class rational {
private:
    T _denominator, _numerator;
    double _value() const {return 1.0*_numerator/_denominator;}
public:
    rational():_denominator(1),_numerator(0){}
    rational(const T& n):_denominator(1),_numerator(n){}
    rational(const T& n, const T& d):_denominator(d),_numerator(n){
        if (_denominator<0){
            _denominator=-_denominator;
            _numerator=-_numerator;
        }
        fraction_reduction(_numerator,_denominator);
    }
    rational(const rational<T> & r):_denominator(r._denominator),_numerator(r._numerator){}
    
    const T denominator() const {return _denominator;}
    const T numerator() const {return _numerator;}
    
    rational<T> &operator=(const rational<T> &r){
        if (this != &r){
            this->_numerator=r._numerator;
            this->_denominator=r._denominator;
        }
    }
    
    rational<T> &operator+=(const rational<T> &r){
        T g=gcd(std::abs(_denominator),std::abs(r._denominator)), s=r._denominator/g;
        this->_denominator*=s;
        this->_numerator*=s;
        this->_numerator+=r._numerator*this->_denominator/g/s;
        fraction_reduction(_numerator,_denominator);
    }
    
    rational<T> &operator-=(const rational<T> &r){
        T g=gcd(_denominator,r._denominator), s=r._denominator/g;
        this->_denominator*=s;
        this->_numerator*=s;
        this->_numerator-=r._numerator*this->_denominator/g/s;
        fraction_reduction(_numerator,_denominator);
    }
    
    rational<T> &operator*=(const rational<T> &r){
        T g1=gcd(std::abs(_denominator),std::abs(r._numerator));
        T g2=gcd(std::abs(r._denominator),std::abs(_numerator));
        
        this->_denominator=(this->_denominator*r._denominator/g1/g2);
        this->_numerator=(this->_numerator*r._numerator/g1/g2);
    }
    
    rational<T> &operator/=(const rational<T> &r){
        return ((*this)*=r.inv());
    }
    
    rational<T> inv() const{
        rational temp;
        temp._denominator=_numerator;
        temp._numerator=_denominator;
        if (temp._denominator<0){
            temp._denominator*=(-1);
            temp._numerator*=(-1);
        }
        return temp;
    }
    inline friend const rational operator+(const rational &M) { return M; }
    
    inline friend const rational operator+(const rational &M, const rational &N) { return rational(M) += N; }
    
    inline friend const rational operator-(const rational &M, const rational &N) { return rational(M) -= N; }
    
    inline friend const rational operator*(const rational &M, const rational &N) { return rational(M) *= N; }
    
    inline friend const rational operator/(const rational &M, const rational &N) { return rational(M) /= N; }
    
    inline friend const rational operator-(const rational &M) { rational temp(M); temp._numerator*=(-1); return temp; }
    
    
    inline friend const bool operator==(const rational<T> & M, const rational<T> & N) {
        return (M._numerator==N._numerator && M._denominator==N._denominator);
    }
    inline friend const bool operator!=(const rational<T> & M, const rational<T> & N) {
        return (M._numerator!=N._numerator || M._denominator!=N._denominator);
    }
    inline friend const bool operator<=(const rational<T> & M, const rational<T> & N) {
        return (M._value()<=N._value());}
    
    inline friend const bool operator>=(const rational<T> & M, const rational<T> & N) {
        return (M._value()>=N._value());}
    
    inline friend const bool operator<(const rational<T> & M, const rational<T> & N) {
        return (M._value()<N._value());}
    
    inline friend const bool operator>(const rational<T> & M, const rational<T> & N) {
        return (M._value()>N._value());}
    
    inline friend const rational abs(const rational &r){
        rational temp(r); temp._numerator=std::abs(temp._numerator); return temp;
    }
    
    
    friend std::ostream& operator<<(std::ostream& output, const rational<T>& r){
        output <<r.numerator()<<"/"<<r.denominator();
        return output;
    }
};


#endif //PYTPSA_RATIONAL_H
