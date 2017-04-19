//
// Created by Yue Hao on 4/4/17.
//
#include "polymap.h"
#include "mathfunc.h"

CPolyMap::CPolyMap(const int& dim, const int& order):dim(dim),max_order(order) {
    setindexmap();
}
std::vector<int> CPolyMap::decomposite(const int& n){
    std::vector<int> result((unsigned long)dim+1);
    int itemp=n+1;
    for (int i=dim; i>0; i--) {
        int k=i-1;
        while (binomial(k, i)<itemp) {
            k++;
        }
        itemp-=binomial(k-1, i);
        result[dim-i]=k-i;
    }
    for (int i=dim;i>0;i--) {
        result[i]=result[i-1]-result[i];
    }
    return result;
}
void CPolyMap::setindexmap(){
    unsigned long totallength=binomial(max_order+dim, dim);
    map.resize(totallength);
    #pragma omp parallel for schedule(dynamic)
    for (int i=0; i<totallength; i++) {
        map[i]=this->decomposite(i);
    }
    return;
}