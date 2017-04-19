

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

#include "tpsalib.h"









/*void CTPS<T>::ErrMsg(const CTPS<T>::errors& err, const string& fromfunc) const{
    switch (err) {
        case OverFlow:
            cout << "The function "<<fromfunc<<" in class CTPS, INDEX OverFlow";
            break;
        case DivZero:
            cout << "The function "<<fromfunc<<" in class CTPS, Divided by zero";
            break;
        case DiffDim:
            cout << "The function "<<fromfunc<<" in class CTPS, has different dimension";
            break;
        case NegValue:
            cout << "The function "<<fromfunc<<" in class CTPS, has negative parameter";
            break;
        default:
            break;
    }
    exit(-1);
    return;
}*/




template<class T>
CTPS<T>::CTPS(){
    //int dim=CTPS<T>::TPS_Dim;
    //degree=1;
    //terms=(unsigned long)dim+1;
    //map.reserve(binomial(dim+CTPS<T>::Maximum_TPS_Degree, dim));
    this->assign(T(0.0));
    
}

template<class T>
CTPS<T>::CTPS(const T& a){
    this->assign(a);
}

template<class T>
CTPS<T>::CTPS(const CTPS<T> &M){
    int dim=TPS_Dim;
    this->degree=M.degree;
    this->terms=M.terms;
    this->map=M.map;


}

template<class T>
unsigned long CTPS<T>::findindex(const std::vector<int>& indexmap) const{
    int dim=TPS_Dim;
    std::vector<int> sum((unsigned long)dim+1);
    sum[0]=indexmap[0];

    for (int i=1; i<=dim; i++) {
        sum[i]=sum[i-1]-indexmap[i];
    }
    unsigned long result=0;
    for (int i=dim; i>0; i--) {
        if (sum[dim-i]==0) {
            break;
        }
        result+=binomial(sum[dim-i]-1+i, i);
    }
    return result;
}

template<class T>
std::vector<int> CTPS<T>::findpower(const unsigned long &n) const {
    return this->polymap.getindexmap(n);
}

template<class T>
void CTPS<T>::redegree(const int& degree){
    this->degree=degree;
    if (degree>CTPS<T>::Maximum_TPS_Degree) this->degree=CTPS<T>::Maximum_TPS_Degree;
    terms=binomial(TPS_Dim+degree, degree);
    map.resize(terms);
}

template <class T>
void CTPS<T>::assign(const T& a, const int& n_var){

    if (n_var<=this->TPS_Dim && n_var>0) {
        this->degree=1;
        this->terms=(unsigned long)CTPS<T>::TPS_Dim+1;
        map.clear();
        //map.reserve(binomial(CTPS<T>::TPS_Dim+CTPS<T>::Maximum_TPS_Degree, CTPS<T>::TPS_Dim));
        map.assign(this->terms, T(0.0));
        this->map[n_var]=T(1.0);
        this->map[0]=a;
    }
    else throw std::runtime_error(std::string("Num of var out of range in CTPS"));
    
}

template <class T>
void CTPS<T>::assign(const T& a){
    this->degree=0;
    this->terms=1;
    map.clear();
    map.reserve(binomial(CTPS<T>::TPS_Dim+CTPS<T>::Maximum_TPS_Degree, CTPS<T>::TPS_Dim));
    map.assign(this->terms, T(0.0));
    this->map[0]=a;
}

template<class T>
const T CTPS<T>::element(const unsigned long & ind) const{
    if (ind<0 || ind >=terms) {
        exit(0);
    }
    return map[ind];
}

template<class T>
const T CTPS<T>::element(const std::vector<int>& ind) const{
    unsigned long result=this->findindex(ind);
    return map[result];
}

template <class T>
T CTPS<T>::evaluate(const std::vector<T>& inivalue) const{
    if (inivalue.size()!=this->TPS_Dim) {
        throw std::runtime_error(std::string("Inconsistent dimension to evaluate CTPS"));
    }
    T sum=map[0];
    #pragma omp parallel for
    for (int i=1; i<this->terms;i++){
        std::vector<int> temp=this->polymap.getindexmap(i);
        T product=T(1.0);
        for (int j=0; j<this->TPS_Dim; j++) {
            product=product*std::pow(inivalue[j], temp[j+1]);
        }
        sum+=product*map[i];
    }
    return sum;
}
template <class T>
CTPS<T> CTPS<T>::derivative(const int& ndim, const int &order) const {
    if (order<=0) return CTPS(*this);
    if (ndim<=this->TPS_Dim && ndim>0) {
        CTPS<T> temp((*this)*T(0.0));
        int new_max_order=0;
        for (int i=1;i<terms;i++){
            std::vector<int> indexlist=this->polymap.getindexmap(i);
            if (indexlist[ndim]>=order){
                int thisdim=indexlist[ndim];
                indexlist[ndim]-=order;
                indexlist[0]-=order;
                if (new_max_order<indexlist[0]) new_max_order=indexlist[0];
                unsigned long new_i=findindex(indexlist);
                temp.map[new_i]=this->map[i]*double(binomial(thisdim, order));
            }
        }
        temp.redegree(new_max_order);
        return temp;
        
    }
    else throw std::runtime_error(std::string("Inconsistent dimension to take derivative"));
    return CTPS();
}

template<class T>
std::ostream& CTPS<T>::print_by_order(std::ostream& output) const {
    int current_order=0;
    int acc_return=0;
    int just_endl=0;
    output<<"Order 0:"<<std::endl;
    for (int i=0;i<this->terms;i++){
        if (std::abs(this->map[i])==0 && i > 0) continue;
        std::vector<int> temp=this->polymap.getindexmap(i);
        if (temp[0] > current_order){
            if (just_endl==0) output<<std::endl;
            output<<std::endl<<"Order "<<temp[0]<<":"<<std::endl;
            current_order=temp[0];
            acc_return=0;
        }
        output<<"(";
        for (int j=1;j<temp.size();j++) {
            output<<temp[j];
            if (j<temp.size()-1) output<<',';
        }
        output<<"): ";
        output.width(10);
        output<<std::left<<this->map[i];
        acc_return++;

        if (acc_return % 3 == 0) {output << std::endl; just_endl=1;}
        else {output << '\t'; just_endl=0;}


    }
    if (just_endl==0) output<<std::endl;
    return output;
}

template<class T>
CTPS<T>& CTPS<T>::operator=(const CTPS<T> &M){
    if (this != &M){
        this->degree=M.degree;
        this->terms=M.terms;
        this->map=M.map;
    }
    return *this;
}

template<class T>
CTPS<T>& CTPS<T>::operator+=(const CTPS<T> &M){

    if (this->degree<M.degree) {
        this->redegree(M.degree);
    }
    //#pragma omp parallel for schedule(static)
    for (int i=0; i<M.terms; i++) {
        this->map[i]+=M.map[i];
    }
    return *this;

}

template<class T>
CTPS<T>& CTPS<T>::operator-=(const CTPS<T>&M){

    if (this->degree<M.degree) {
        this->redegree(M.degree);
    }
    //#pragma omp parallel for schedule(static)
    for (int i=0; i<M.terms; i++) {
        this->map[i]-=M.map[i];
    }
    return *this;
}

template<class T>
CTPS<T>& CTPS<T>::operator*=(const CTPS<T>& M){

    if (M.get_degree()==0){

        //#pragma omp parallel for schedule(static,4096)
        for (int i=0; i<this->terms; i++){
            this->map[i]*=M.map[0];
        }
        return *this;
    }

    CTPS<T>temp(*this);
    this->map.clear();
    (*this).redegree(std::min(CTPS<T>::Maximum_TPS_Degree, this->degree+M.degree));
    //#pragma omp parallel for schedule(dynamic)
    for (int i=0; i<temp.map.size(); i++ ) {
        if (std::abs(temp.map[i])==0) {
            continue;
        }
        std::vector<int> vthis=polymap.getindexmap(i);
        //#pragma omp parallel for schedule(dynamic)
        for (int j=0; j<M.map.size();j++) {
            if (std::abs(M.map[j])==0) {
                continue;
            }
            std::vector<int> vm=polymap.getindexmap(j);
            if (vthis[0]+vm[0]>CTPS<T>::Maximum_TPS_Degree) break;
            std::vector<int> indexmap(vthis.size());
            for (int k=0; k<vthis.size(); k++) {
                indexmap[k] = vthis[k] + vm[k];
            }
            unsigned long target_ind=this->findindex(indexmap);
            //#pragma omp atomic
            this->map[target_ind]+=temp.map[i]*M.map[j];
        }
    }
    return *this;
}

template<class T>
CTPS<T>& CTPS<T>::operator/=(const CTPS<T>& M){

    if (abs(M.cst())==0) {
        throw std::runtime_error("Divide by zero, in CTPS");
    }
    if (M.get_degree()==0){
        #pragma omp parallel for schedule(static)
        for (int i=0; i<this->terms; i++){
            this->map[i]/=M.map[0];
        }
        return *this;
    }
    return ((*this)*=inv(M));


}



