

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
CTPS<T>::CTPS(const T& a, const int& n){
    this->assign(a,n);
}

template<class T>
CTPS<T>::CTPS(const CTPS<T> &M){
    int dim=TPS_Dim;
    this->degree=M.degree;
    this->terms=M.terms;
    this->map=M.map;
    //this->index=M.index;

}
template<class T>
CTPS<T>::CTPS(const std::vector<T>& input_map){
    this->set_map();

}

template<class T>
unsigned long CTPS<T>::find_index(const std::vector<int>& indexmap) const{
    int dim=TPS_Dim;
    if (indexmap.size()==dim){
        std::vector<int> newindexmap((unsigned long)dim+1);
        int sum=0;
        for (int i=1;i<=dim;i++){
            sum+=indexmap[i-1];
            newindexmap[i]=indexmap[i-1];
        }
        newindexmap[0]=sum;
        return this->find_index(newindexmap);

    }
    if (indexmap.size()!= (dim+1)) throw std::runtime_error(std::string("Index map does not have correction length"));
    std::vector<int> sum((unsigned long)dim+1);
    sum[0]=indexmap[0];

    for (int i=1; i<=dim; i++) {
        if (indexmap[i]<0) throw std::runtime_error(std::string("The index map has invalid component"));
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
std::vector<int> CTPS<T>::find_power(const unsigned long &n) const {
    if (n < this->Maximum_Terms)  return this->polymap.getindexmap(n);
    else
        throw std::out_of_range("The index is out of range");
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
        map.reserve(binomial(CTPS<T>::TPS_Dim+CTPS<T>::Maximum_TPS_Degree, CTPS<T>::TPS_Dim));
        map.assign(this->terms, T(0.0));
        this->map[n_var]=T(1.0);
        this->map[0]=a;
        //this->index.insert(0);this->index.insert(n_var);
    }
    else throw std::runtime_error(std::string("Num of var out of range in CTPS"));
    
}

template <class T>
void CTPS<T>::assign(const T& a){
    this->degree=0;
    this->terms=1;
    this->map.clear();
    this->map.push_back(a);
    //map.reserve(binomial(CTPS<T>::TPS_Dim+CTPS<T>::Maximum_TPS_Degree, CTPS<T>::TPS_Dim));
    
    //this->index.insert(0);
}

template<class T>
const T CTPS<T>::element(const unsigned long & ind) const{
    if (ind<0 || ind >=terms) {
        //throw std::runtime_error(std::string("Element index out of range in CTPS"));
        throw std::out_of_range(std::string("Element index out of range in CTPS"));
    }
    return map[ind];
}

template<class T>
void CTPS<T>::set_map(const std::vector<T>& input_map){
    unsigned ord;
    unsigned long term_of_ord;
    for (ord=0; ord < CTPS<T>::Maximum_TPS_Degree+1; ord ++){
        term_of_ord=binomial(CTPS<T>::TPS_Dim+ord, CTPS<T>::TPS_Dim);

        if (input_map.size()>term_of_ord){
            continue;
        }
    }
    //std::cout<<ord<<'\t'<<term_of_ord<<'\t'<<input_map.size()<<std::endl;
    this->degree=ord;
    this->terms=term_of_ord;
    this->map=input_map;
    this->map.resize(this->terms, T(0.0));

}

template<class T>
const T CTPS<T>::element(const std::vector<int>& ind) const{
    unsigned long result=this->find_index(ind);
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
CTPS<T> CTPS<T>::evaluate(const std::vector<CTPS<T> >& inivalue) const{
    if (inivalue.size()!=this->TPS_Dim) {
        throw std::runtime_error(std::string("Inconsistent dimension to evaluate CTPS"));
    }
    CTPS<T> sum(map[0]);
    #pragma omp parallel for
    for (int i=1; i<this->terms;i++){
        std::vector<int> temp=this->polymap.getindexmap(i);
        CTPS<T> product(T(1.0));
        for (int j=0; j<this->TPS_Dim; j++) {
            for (int k=0; k<temp[j+1]; k++)
            product*=inivalue[j];
        }
        sum+=product*map[i];
    }
    return sum;
}

template<class T>
CTPS<T> CTPS<T>::linear() const{
    CTPS<T> result(T(0.0));
    result.degree=1;
    result.terms=(unsigned long)TPS_Dim+1;
    result.map.resize(result.terms);
    for (unsigned long i=0;i<std::min(result.terms,this->terms);i++){
        //std::cout<<i<<std::endl;
        result.map[i]=this->map[i];
    }
    return result;

}

template<class T>
bool CTPS<T>::is_zero() const{

    for (unsigned long i = 0; i < this->get_terms(); i++)
        if (std::abs(this->map[i]) > inf_epsilon) return false;
    return true;

}

template<class T>
bool CTPS<T>::is_equal(const CTPS<T> &M) const {
    unsigned long this_terms = this->get_terms();
    unsigned long m_terms = M.get_terms();
    unsigned long terms = this_terms;

    if (this_terms > m_terms) terms = m_terms;
    for (unsigned long i = 0; i < terms; i++) {
        if (std::abs(this->map[i] - M.map[i]) < inf_epsilon ||
            std::abs(this->map[i] - M.map[i]) < inf_epsilon * std::abs(this->map[i])) {}
        else { return false; }
    }

    if (this_terms > m_terms) {
        for (unsigned long i = terms; i < this_terms; i++) {
            if (std::abs(this->map[i]) > inf_epsilon) return false;
        }
    }
    if (this_terms < m_terms) {
        for (unsigned long i = terms; i < m_terms; i++) {
            if (std::abs(M.map[i]) > inf_epsilon) return false;
        }
    }


}

template <class T>
CTPS<T> CTPS<T>::derivative(const int& ndim, const int &order) const {
    if (order<=0) return CTPS(*this);
    if (this->get_degree()<order) return CTPS(T(0.0));
    if (ndim<=this->TPS_Dim && ndim>0) {
        CTPS<T> temp((*this)*T(0.0));
        int new_max_order=0;
        for (int i=1;i<terms;i++){
            std::vector<int> indexlist=this->polymap.getindexmap(i);
            if (indexlist[ndim]>=order){
                int thisdim=indexlist[ndim];
                indexlist[ndim]-=order;
                indexlist[0]-=order;
                //if (new_max_order<indexlist[0]) new_max_order=indexlist[0];
                unsigned long new_i=find_index(indexlist);
                temp.map[new_i]=this->map[i]*double(binomial(thisdim, order));
            }
        }
        temp.redegree(this->get_degree()-order);
        return temp;
        
    }
    else throw std::runtime_error(std::string("Inconsistent dimension to take derivative"));
    return CTPS();
}

template <class T>
CTPS<T> CTPS<T>::integrate(const int& ndim, const T &a0) const {
    if (ndim<=this->TPS_Dim && ndim>0) {
        CTPS<T> temp((*this)*T(a0));
        int new_max_order=this->get_degree()+1;
        temp.redegree(new_max_order);
        for (int i=0;i<terms;i++){
            std::vector<int> indexlist=this->polymap.getindexmap(i);
            if (indexlist[ndim]<new_max_order){
                int thisdim=indexlist[ndim];
                indexlist[ndim]+=1;
                indexlist[0]+=1;
                unsigned long new_i=find_index(indexlist);
                temp.map[new_i]=this->map[i]/(thisdim+1.0);
            }
        }

        return temp;

    }
    else throw std::runtime_error(std::string("Inconsistent dimension to take integration"));
    return CTPS();
}

template <class T>
CTPS<T> CTPS<T>::conjugate(const int &mode) const{
    CTPS<T> temp(*this);
    if (this->TPS_Dim % 2 !=0 ) return temp;
    std::vector<int> nid, ncid;
    nid.resize(this->TPS_Dim / 2, 0);
    ncid.resize(this->TPS_Dim / 2, 0);
    if (mode==1) { //z1, z1*, z2, z2* ...
        for (int i = 1; i<=this->TPS_Dim / 2; i++){
            nid[i-1]=i*2-1;
            ncid[i-1]=i*2;
        }
    }
    if (mode==2) { //z1, z2, z3..., z1*, z2*, z2*...
        for (int i = 1; i<=this->TPS_Dim / 2; i++){
            nid[i-1]=i;
            ncid[i-1]=i+this->TPS_Dim / 2;
        }
    }
    for (unsigned long j = 1; j<this->terms; j++){
            std::vector<int> vthis=polymap.getindexmap(j);
            for (int i = 0; i<nid.size(); i++){
                int temp=vthis[nid[i]];
                vthis[nid[i]]=vthis[ncid[i]];
                vthis[ncid[i]]=temp;
            }
            unsigned long newj=this->find_index(vthis);
            temp.map[newj]=(this->map[j]);
        }
        return temp;
}


template <>
CTPS<std::complex<double>> CTPS<std::complex<double>>::conjugate(const int &mode) const{
    CTPS<std::complex<double> > temp(*this);
    if (this->TPS_Dim % 2 !=0 ) return temp;
    std::vector<int> nid, ncid;
    nid.resize(this->TPS_Dim / 2, 0);
    ncid.resize(this->TPS_Dim / 2, 0);
    if (mode==1) { //z1, z1*, z2, z2* ...
        for (int i = 1; i<=this->TPS_Dim / 2; i++){
            nid[i-1]=i*2-1;
            ncid[i-1]=i*2;
        }
    }
    if (mode==2) { //z1, z1*, z2, z2* ...
        for (int i = 1; i<=this->TPS_Dim / 2; i++){
            nid[i-1]=i;
            ncid[i-1]=i+this->TPS_Dim / 2;
        }
    }
    for (unsigned long j = 1; j<this->terms; j++){
            std::vector<int> vthis=polymap.getindexmap(j);
            for (int i = 0; i<nid.size(); i++){
                int temp=vthis[nid[i]];
                vthis[nid[i]]=vthis[ncid[i]];
                vthis[ncid[i]]=temp;
            }
            unsigned long newj=this->find_index(vthis);

            temp.map[newj]=std::conj(this->map[j]);
        }
    return temp;
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
        //this->index=M.index;
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
    if (this==&M) {
        CTPS<T> mcopy(M);
        return (*this)*=mcopy;
    }
    //std::vector<T> mmap=M.map;
    //unsigned long msize=M.map.size();
    int newdegree=this->degree+M.degree;


    CTPS<T>temp(*this);
    this->map.clear();
    (*this).redegree(std::min(CTPS<T>::Maximum_TPS_Degree, newdegree));
    //#pragma omp parallel for schedule(dynamic)
    for (int i=0; i<temp.map.size(); i++ ) {
        if (std::abs(temp.map[i])==0) {
            continue;
        }
        std::vector<int> vthis=polymap.getindexmap(i);
        //#pragma omp parallel for schedule(dynamic)
        //for (int j=0; j<msize();j++) {
        unsigned long j_max=std::min(M.map.size(),binomial(this->TPS_Dim+this->Maximum_TPS_Degree-vthis[0],this->TPS_Dim));
        //#pragma omp parallel for schedule(dynamic)
        for (int j=0; j< j_max;j++) {
            if (std::abs(M.map[j])==0) {
                continue;
            }
            std::vector<int> vm=polymap.getindexmap(j);
            //if (vthis[0]+vm[0]>CTPS<T>::Maximum_TPS_Degree) break;
            std::vector<int> indexmap(vthis.size());
            for (int k=0; k<vthis.size(); k++) {
                indexmap[k] = vthis[k] + vm[k];
            }
            unsigned long target_ind=this->find_index(indexmap);
            //#pragma omp atomic
            this->map[target_ind]+=temp.map[i]*M.map[j];
        }
    }
    return *this;
}

template<class T>
CTPS<T>& CTPS<T>::operator/=(const CTPS<T>& M){

    if (std::abs(M.cst())==0) {
        throw std::runtime_error("Divide by zero, in CTPS");
    }
    if (this==&M) {
        (*this)=CTPS<T>(T(1.0));
        return *this;
    }
    if (M.get_degree()==0){
        //#pragma omp parallel for schedule(static)
        for (int i=0; i<this->terms; i++){
            this->map[i]/=M.map[0];
        }
        return *this;
    }
    return ((*this)*=inv(M));


}



