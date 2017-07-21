%module tpsalib
%{
  #define SWIG_FILE_WITH_INIT
  #include "../src/tpsalib.h"
  #include <sstream>
  #include <string>
  #include <complex>
%}

%include <std_vector.i>
%include <std_complex.i>
%include <stl.i>
%include <std_string.i>


%template(DoubleVector) std::vector<double>;
%template(IntVector) std::vector<int>;



%include "../src/tpsalib.h"


%extend CTPS{

  CTPS<T> __add__(const CTPS<T>& a){
    return CTPS<T>(*$self)+=a;
  }

  CTPS<T> __add__(const T& a){
    return CTPS<T>(*$self)+=T(a);
  }
  CTPS<T> __radd__(const T& a){
    return CTPS<T>(*$self)+=T(a);
  }

  CTPS<T> __sub__(const CTPS<T>& a){
    return CTPS<T>(*$self)-=a;
  }
  CTPS<T> __sub__(const T& a){
    return CTPS<T>(*$self)-=T(a);
  }
  CTPS<T> __rsub__(const T& a){
    return CTPS<T>(*$self)*T(-1.0)+T(a);
  }

  CTPS<T> __mul__(const CTPS<T>& a){
    return CTPS<T>(*$self)*=a;
  }
  CTPS<T> __mul__(const T& a){
    return CTPS<T>(*$self)*=T(a);
  }
  CTPS<T> __rmul__(const T& a){
    return CTPS<T>(*$self)*=T(a);
  }

  CTPS<T> __div__(const CTPS<T>& a){
    return CTPS<T>(*$self)*=inv(a);
  }
  CTPS<T> __div__(const T& a){
    return CTPS<T>(*$self)/=T(a);
  }
  CTPS<T> __rdiv__(const T& a){
    return inv(*$self)*T(a);
  }
  CTPS<T> __neg__(){
    return CTPS<T>(*$self)*T(-1.0);
  }
  CTPS<T> __pos__(){
    return CTPS<T>(*$self);
  }
  
  std::string  __repr__(){
    std::ostringstream out;
    $self->print_by_order(out);
    std::string output=out.str();
      return output;
  }
 

}


%ignore operator<<;
%ignore operator=;

/*%rename (__add__) operator+;
%rename (__sub__) operator-;
%rename (__div__) operator/;
%rename (__mul__) operator*;

%rename (__neg__) operator-();
%rename (__pos__) operator+();
%rename (__eq__) operator==;
%rename (__ne__) operator!=;*/

%template(DCTPS) CTPS<double>;
%template(CCTPS) CTPS<std::complex<double> >;

