//
// Created by Yue Hao on 4/5/17.
//

#include "../tpsalib.h"
#include <complex>
#include <iostream>
int main(){
    cctps::Initialize(6,6);
    cctps x,y,z;
    x.assign(1.0+1i,1);
    y.assign(0-1i,2);
    z=inv(x+y);
    std::cout<<z;
}