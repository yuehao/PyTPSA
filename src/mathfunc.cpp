//
// Created by Yue Hao on 4/4/17.
//


unsigned long factorial(const int & n) {
    if (n<0) {
        return 0;
    }
    switch (n) {
        case 0:
            return 1;
        case 1:
            return 1;
        case 2:
            return 2;
        case 3:
            return 6;
        case 4:
            return 24;
        case 5:
            return 120;
        case 6:
            return 720;
        case 7:
            return 5040;
        case 8:
            return 40320;
        case 9:
            return 362880;
        case 10:
            return 3628800;
        case 11:
            return 39916800;
        case 12:
            return 479001600;
        
        
        default:
            return n*factorial(n-1);
    }
}
unsigned long doublefactorial(const int& n){
    if (n<0) {
        return 0;
    }
    switch (n) {
        case 0:
            return 1;
        case 1:
            return 1;
        case 2:
            return 2;
        case 3:
            return 3;
        case 4:
            return 8;
        case 5:
            return 15;
        case 6:
            return 48;
        case 7:
            return 105;
        case 8:
            return 384;
        case 9:
            return 945;
        case 10:
            return 3840;
        case 11:
            return 10395;
        case 12:
            return 46080;
        case 13:
            return 135135;
        case 14:
            return 645120;
        case 15:
            return 2027025;
        case 16:
            return 10321920;
        case 17:
            return 34459425;
        case 18:
            return 185794560;
        
        default:
            return n*doublefactorial(n-2);
    }
    
}
unsigned long binomial(const int& n, const int& m) {
    if (n<=0 || m>n || m<0) {
        return 0;
    }
    int ml;
    if (m>n/2) ml=n-m;
    else ml=m;
    switch (ml) {
        case 0:
            return 1;
        case 1:
            return (unsigned long)n;
        case 2:
            return (unsigned long)n*(n-1)/2;
        
        case 3:
            return (unsigned long)n*(n-1)*(n-2)/6;
        
        
        default:
            break;
    }
    if (n<=12) {
        return factorial(n)/factorial(m)/factorial(n-m);
    }
    else {
        return (n*binomial(n-1,ml-1)/ml);
    }
}