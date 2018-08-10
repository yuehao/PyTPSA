//
// Created by Yue Hao on 4/4/17.
//

#include "../tpsalib.h"
#include <iostream>
#include <vector>
#include <ctime>

int main(){
    dctps::Initialize(8,10);
    dctps x,y;
    x.assign(1.0,1);
    x*=x;
    x*=x;
    x/=x;
    std::cout<<x;
    
    /*for (int i=0;i<1;i++){
      std::clock_t start;
      start=std::clock();
      dctps x,y;
      std::vector<dctps> z;
      x.assign(1.0,1);
      y.assign(1.0,2);
      int zsize=100000;
      z.resize(zsize);
      for (int j=0;j<zsize;j++) z[i]=x;
      for (int k=0;k<10;k++)
	for (int j=0;j<zsize;j++) z[i]=(z[i]-1.0)*(z[i]+1.0);
      
      double duration=(std::clock()-start)/(double)CLOCKS_PER_SEC;
      std::cout<<duration<<std::endl;
      //z=inv(x+y);
      //x=z.derivative(1,1);
    //std::cout<<z;
    }*/
}
