//
// Created by Yue Hao on 4/4/17.
//

#include "../tpsalib.h"
#include <iostream>

int main(){
    dctps::Initialize(6,6);
    dctps x,y,z;
    x.assign(1.0,1);
    y.assign(0.0,2);
    z=inv(x+y);
    std::cout<<z;
    
}