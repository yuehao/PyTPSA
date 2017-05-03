//
// Created by Yue Hao on 4/4/17.
//

#include "../tpsalib.h"
#include <iostream>

int main(){
    dctps::Initialize(6,6);
    for (int i=0;i<100;i++){
    dctps x,y,z;
    x.assign(1.0,1);
    y.assign(0.0,2);
    x=1+x*(y+1);
    
    z=inv(x+y);
    x=z.derivative(1,1);
    //std::cout<<z;
    }
}
