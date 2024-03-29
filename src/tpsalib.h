/*
 *  tpsalib.h
 *
 *  Created by Yue Hao on 4/8/11.
 *
 */


#ifndef TPS
#define TPS
#include <complex>
#include <vector>
#include <iostream>
#include <sstream>
#include <string>
#include <set>
#include "polymap.h"
#include "mathfunc.h"
#include "rational.h"
//using namespace std;
const int MAX_TPS_ORDERS=6;
const double inf_epsilon=1e-16;



template<class T> // T is the datatype for the matrix
class CTPS {
private:
    int degree;
    unsigned long terms;//Total terms
    std::vector<T> map;
    std::set<unsigned long> indices;
    void redegree(const int &);

public:
    static int Maximum_TPS_Degree;
    static int TPS_Dim;
    static int Maximum_Terms;
    static CPolyMap polymap;
    
    static void Initialize(const int &dim, const int &max_order = MAX_TPS_ORDERS) {
        CTPS<T>::TPS_Dim = dim;
        CTPS<T>::Maximum_TPS_Degree = max_order;
        CTPS<T>::Maximum_Terms = binomial(dim+max_order, dim);
        CTPS<T>::polymap = CPolyMap(dim, max_order);
    }
    static int Get_Max_Degree(){
        return Maximum_TPS_Degree;
    }
    static int Get_TPS_Dim(){
        return TPS_Dim;
    }
    static int Get_Max_Terms(){
        return Maximum_Terms;
    }
    static std::vector<int> Get_Power_Index(const int & i){
        return polymap.getindexmap(i);
    }
    
    CTPS();
    
    CTPS(const T &a);  //Constant
    
    CTPS(const T &a, const int &);  //Variable
    
    CTPS(const CTPS<T> &);

    CTPS(const std::vector<T>& input_map);
    
    ~CTPS() {}
    
    void assign(const T &); //A constant
    void assign(const T &, const int &);// A Variable
    
    unsigned long find_index(const std::vector<int> &) const;
    
    std::vector<int> find_power(const unsigned long &n) const;
    
    inline const int get_dim() const { return TPS_Dim; }
    
    inline const int get_degree() const { return degree; }
    
    inline const unsigned long get_terms() const { return terms; }
    
    const T element(const unsigned long &ind) const;
    
    const T element(const std::vector<int> &indmap) const;

    const std::vector<T> get_map() const {return this->map;}
    void set_map(const std::vector<T>& input_map);

    void machine_precision_roundoff();

    std::ostream&  print_by_order(std::ostream& output) const;
    
    //operator double() const;
    CTPS<T> &operator=(const CTPS<T> &);
    
    CTPS<T> &operator+=(const CTPS<T> &);
    
    CTPS<T> &operator-=(const CTPS<T> &);
    
    CTPS<T> &operator*=(const CTPS<T> &);
    
    CTPS<T> &operator/=(const CTPS<T> &);
    
    CTPS<T> &add_to(const CTPS<T> &M) {(*this) += M; return *this;}
    
    CTPS<T> &sub_to(const CTPS<T> &M) { (*this) -= M; return *this;}
    
    CTPS<T> &mul_to(const CTPS<T> &M) { (*this) *= M; return *this;}
    
    CTPS<T> &div_to(const CTPS<T> &M) { (*this) /= M; return *this;}

    bool is_zero() const;
    bool is_equal(const CTPS<T> &M) const;
    
    T evaluate(const std::vector<T> &value) const;
    CTPS<T> evaluate(const std::vector<CTPS<T> > &) const;
    
    CTPS<T> derivative(const int &ndim, const int &order = 1) const;
    CTPS<T> integrate(const int &ndim, const T& ) const;
    CTPS<T> conjugate(const int &mode=1) const;

    
    inline const T cst() const { return map[0]; }
    CTPS<T> linear() const;

    template<class U> CTPS<U> convert_to() const{
        CTPS<U> temp;
        temp.degree=this->degree;
        temp.indices=this->indices;
        temp.terms=this->terms;
        temp.map.resize(temp.terms);
        for (int i=0;i<this->terms; i++){
            temp.map[i]=U(this->map[i]);
        }
    }

    
    
    friend CTPS inv(const CTPS & M) {
        CTPS temp(M), sum, term_by_order;
        if (std::abs(M.cst()) == 0) {
            throw std::runtime_error("Divide by zero, in CTPS");
        }
        temp = temp - M.cst();
        term_by_order = T(1) / M.cst();
        sum += term_by_order;
        for (int i = 1; i<=M.Maximum_TPS_Degree; i++) {
            term_by_order *= (-temp / M.cst());
            sum += term_by_order;
        }
        return sum;
    }
    
    friend CTPS exp(const CTPS & M) {
        CTPS temp(M), sum, term_by_order;
        //T a0=M.cst();
        temp = temp - M.cst();
        term_by_order += T(1);
        sum += term_by_order;
        for (int i = 1; i <= M.Maximum_TPS_Degree; i++) {
            double index = 1.0 / factorial(i);
            term_by_order *= temp;
            sum += (term_by_order * T(index));
        }
        sum *= exp(M.cst());
        return sum;
    }
    
    friend CTPS log(const CTPS & M) {
        CTPS temp(M), sum, term_by_order;
        temp = temp - M.cst();
        term_by_order = temp / M.cst();
        sum += term_by_order;
        for (int i = 2; i <= M.Maximum_TPS_Degree; i++) {
            //T index=pow(-1.0, i+1)/i/pow(a0, i);
            term_by_order *= (-temp / M.cst());
            sum += (term_by_order / T(i));
        }
        sum += log(M.cst());
        return sum;
    }
    
    friend CTPS sqrt(const CTPS & M) {
        CTPS temp(M), sum, term_by_order;
        T a0=sqrt(M.cst());
        //T value=sqrt(a0);
        temp = temp - M.cst();
        term_by_order = temp / a0;
        sum += term_by_order / T(2);
        for (int i = 2; i <= M.Maximum_TPS_Degree; i++) {
            //T index=pow(-1.0, i+1)*doublefactorial(2*i-3)/pow(a0, i-0.5)/doublefactorial(2*i);
            T index = T(1.0 * doublefactorial(2 * i - 3) / doublefactorial(2 * i));
            term_by_order = (-temp) * term_by_order / M.cst();
            sum += (term_by_order * index);
        }
        sum += a0;
        return sum;
    }
    
    friend CTPS pow(const CTPS & M, const double & b) {
        if (b == 1.0) return M;
        if (b == 0.0) {
            CTPS temp;
            return temp += T(1);
        }
        CTPS temp(M), sum, term_by_order;
        double index = b;

        if (M.cst() == T(0.0)){
            if (std::fmod(b, 1.0) == 0.0 && b > 1){
                sum = M;
                for (int i=1; i < b; i++ ) sum = sum * M;
                return sum;
            }
            else{
                std::runtime_error("Divide by zero, in CTPS::pow");
            }
        }
    
        temp = temp - M.cst();
        term_by_order += T(1);
        T factor = pow(M.cst(), b);
        sum += factor;
    
        for (int i = 1; i <= M.Maximum_TPS_Degree; i++) {
            factor = factor / M.cst() * T(index / i);
            index--;
            term_by_order = term_by_order * temp;
            sum += (term_by_order * factor);
            if (index == 0.0) break;
        }
        return sum;
    }
    
    friend CTPS sin(const CTPS & M) {
        CTPS temp(M), sum, term_by_order;
        T a0 = M.cst(), sin_a0 = sin(a0), cos_a0 = cos(a0);
        temp = temp - a0;
        term_by_order += T(1);
        for (int i = 1; i <= M.Maximum_TPS_Degree; i++) {
            T index;
            if (i % 2 == 1) {
                index = cos_a0 * T(pow(-1.0, (i - 1) / 2) / factorial(i));
            } else {
                index = sin_a0 * T(pow(-1.0, i / 2) / factorial(i));
            }
            term_by_order = term_by_order * temp;
            sum += (term_by_order * index);
        }
        sum += sin_a0;
        return sum;
    }

    friend CTPS arcsin(const CTPS & M) {
        /*  Provided by Derong, using fix point trick in M.Berz's book
            y0+y=arcsin(x0+x) sin(y0)=x0
            sin(y0)cos(y)+cos(y0)sin(y)=x0+x
            sin(y)=[x+x0*(1-cos(y))]/cos(y0)
            y=y-sin(y)+[x+x0*(1-cos(y))]/cos(y0)
        */
        CTPS temp(M), result, temp1;
        T a0 = M.cst(), arcsin_a0 = asin(a0), cos_y0=sqrt(T(1.0)-a0*a0);
        temp = temp - a0;
        temp1=temp;
        for (int i=0;i<M.Get_Max_Degree()+1;i++){
            temp1=temp1-sin(temp1)+(temp+a0*(T(1.0)-cos(temp1)))/cos_y0;
        }
        return temp1+arcsin_a0;


    }


    friend CTPS cos(const CTPS & M) {
        CTPS temp(M), sum, term_by_order;
        T a0 = M.cst(), sin_a0 = sin(a0), cos_a0 = cos(a0);
        temp = temp - a0;
        term_by_order += T(1);
        for (int i = 1; i <= M.Maximum_TPS_Degree; i++) {
            T index;
            if (i % 2 == 1) {
                index = sin_a0 * T(pow(-1.0, (i + 1) / 2) / factorial(i));
            } else {
                index = cos_a0 * T(pow(-1.0, i / 2) / factorial(i));
            }
            term_by_order = term_by_order * temp;
            sum += (term_by_order * index);
        }
        sum += cos_a0;
        return sum;
    }

    friend CTPS arccos(const CTPS & M) {
        return T(3.14159265358979323846/2.0)-arcsin(M);
    }

    friend CTPS tan(const CTPS & M) {return sin(M)/cos(M);}


    friend CTPS sinh(const CTPS & M) {
        CTPS temp(M), sum, term_by_order;
        T a0 = M.cst(), sinh_a0 = sinh(a0), cosh_a0 = cosh(a0);
        temp = temp - a0;
        term_by_order += T(1);
        for (int i = 1; i <= M.Maximum_TPS_Degree; i++) {
            T index;
            if (i % 2 == 1) {
                index = cosh_a0 / T(1.0*factorial(i));
            } else {
                index = sinh_a0 / T(1.0*factorial(i));
            }
            term_by_order = term_by_order * temp;
            sum += (term_by_order * index);
        
        }
        sum += sinh_a0;
        return sum;
    }
    
    friend CTPS cosh(const CTPS & M) {
        CTPS temp(M), sum, term_by_order;
        T a0 = M.cst(), sinh_a0 = sinh(a0), cosh_a0 = cosh(a0);
        temp = temp - a0;
        term_by_order += T(1);
        for (int i = 1; i <= M.Maximum_TPS_Degree; i++) {
            T index;
            if (i % 2 == 1) {
                index = sinh_a0 / T(1.0*factorial(i));
            } else {
                index = cosh_a0 / T(1.0*factorial(i));
            }
            term_by_order = term_by_order * temp;
            sum += (term_by_order * index);
        }
        sum += cosh_a0;
        return sum;
    }
    
    inline friend const CTPS operator+(const CTPS &M) { return M; }
    
    inline friend const CTPS operator+(const CTPS &M, const CTPS &N) {
        if (M.get_degree() > N.get_degree())
            return CTPS(M) += N;
        else return CTPS(N) += M;
    }
    
    inline friend const CTPS operator-(const CTPS &M, const CTPS &N) {
        if (M.get_degree() > N.get_degree())
            return CTPS(M) -= N;
        else return CTPS(-N) += M;
    }
    
    inline friend const CTPS operator*(const CTPS &M, const CTPS &N) {
        if (M.get_degree() > N.get_degree())
            return CTPS(M) *= N;
        else return CTPS(N) *= M;
    }
    
    inline friend const CTPS operator/(const CTPS &M, const CTPS &N) { return CTPS(M) /= N; }
    
    inline friend const CTPS operator-(const CTPS &M) { return M * T(-1.0); }
    friend std::ostream& operator<<(std::ostream& output, const CTPS& A){A.print_by_order(output); return output;}
    inline friend const bool operator==(const CTPS & M, const CTPS & N){return M.map==N.map;}
    inline friend const bool operator!=(const CTPS & M, const CTPS & N){return M.map!=N.map;}

    std::string print_to_string() const {std::ostringstream ost; this->print_by_order(ost); return ost.str();}
    
};


using dctps=CTPS<double>;
using cctps=CTPS<std::complex<double> >;
using rctps=CTPS<rational<int> >;
using lrctps=CTPS<rational<long> >;



template<class T>
int CTPS<T>::Maximum_TPS_Degree=-1;
template<class T>
int CTPS<T>::Maximum_Terms=-1;
template<class T>
int CTPS<T>::TPS_Dim=-1;
template<class T>
CPolyMap CTPS<T>::polymap=CPolyMap(1,1);
    
typedef std::complex<double> cplx;

    
    
    /*template<class T>
    inline const CTPS<T> operator+(const CTPS<T> & M) {return M;}
    template<class T>
    inline const CTPS<T> operator-(const CTPS<T> & M) {return CTPS<T>(M)*=(-1.0);}
    template<class T>
    inline const CTPS<T> operator+(const CTPS<T> & M, const CTPS<T> & N) {if (M.get_degree()>N.get_degree()) return CTPS<T>(M)+=N;else return CTPS<T>(N)+=M;}
    template<class T>
    inline const CTPS<T> operator-(const CTPS<T> & M, const CTPS<T> & N) {if (M.get_degree()>N.get_degree()) return CTPS<T>(M)-=N;else return CTPS<T>(-N)+=M;}
    template<class T>
    inline const CTPS<T> operator*(const CTPS<T> & M, const CTPS<T> & N) {if (M.get_degree()>N.get_degree()) return CTPS<T>(M)*=N;else return CTPS<T>(N)*=M;}
    template<class T>
    inline const CTPS<T> operator/(const CTPS<T> & M, const CTPS<T> & N) {return CTPS<T>(M)/=N;}
    //*/




#include "tpsalib.cpp"


#endif
