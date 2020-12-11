cimport cython
cimport libcpp.complex

from libcpp.vector cimport vector
from libcpp.string cimport string

from libcpp cimport bool

cimport tpsa_def

ctypedef double complex cplx


cdef class PyDTPSA:
    cdef tpsa_def.CTPS[double] tps

    def __cinit__(self, value=0.0, int dim_ind=0, PyDTPSA copy_tpsa=None, vector[double] input_map=[]):
        if copy_tpsa is not None:
            self.tps=copy_tpsa.tps
            return

        if tpsa_def.CTPS[double].Get_Max_Degree()<0:
            print("The class is not initialized yet.")
            return

        if len(input_map)>0:
            self.tps.set_map(input_map)
            return

        if isinstance(value,(float,int)):
            if dim_ind>0:
                if dim_ind <=tpsa_def.CTPS[double].Get_TPS_Dim():
                    self.tps.assign(float(value), int(dim_ind))
                else:
                    raise IndexError ("TPSA dimension is out of range.")
            else:
                self.tps.assign(float(value))


    @classmethod
    def initialize(cls, int dim, int max_order):
        tpsa_def.CTPS[double].Initialize(dim, max_order)

    @classmethod
    def get_power_index(cls, int ind):
        return tpsa_def.CTPS[double].Get_Power_Index(ind)

    @classmethod
    def get_max_degree(cls):
        return tpsa_def.CTPS[double].Get_Max_Degree()

    @classmethod
    def get_max_terms(cls):
        return tpsa_def.CTPS[double].Get_Max_Terms()

    def get_dim(self):
        return self.tps.get_dim()
    def get_degree(self):
        return self.tps.get_degree()
    def get_term(self):
        return self.tps.get_terms()

    def indices(self):
        cdef vector[double] tpsmap=self.tps.get_map()
        return tpsmap


    def find_index(self, vector[int] power_ind):
        cdef unsigned long temp
        temp=self.tps.find_index(power_ind)
        return temp

    def find_power(self, unsigned long n):
        cdef vector[int] temp = self.tps.find_power(n)
        return temp

    #def evaluate(self, vector[double] x):
    #    return self.tps.evaluate(x)

    def evaluate(self, vector[double] values):
        return self.tps.evaluate(values)


    def composition(self, values):
        cdef vector[tpsa_def.CTPS[double]] vect
        cdef PyDTPSA a
        for a in values:
            vect.push_back(a.tps)
        cdef PyDTPSA result
        result=PyDTPSA(0.0)
        result.tps=self.tps.evaluate(vect)
        return result



    def derivative(self, int variable, int order):
        cdef PyDTPSA result
        result=PyDTPSA(0.0)
        result.tps=self.tps.derivative(variable, order)
        return result

    def integrate(self, int variable, double a0):
        cdef PyDTPSA result
        result=PyDTPSA(0.0)
        result.tps=self.tps.integrate(variable, a0)
        return result

    def conjugate(self, int mode=1):
        cdef PyDTPSA result
        result=PyDTPSA(0.0)
        result.tps=self.tps.conjugate(mode)
        return result

    def element(self, *l):
        cdef vector[int] indtuple
        cdef int ind
        if len(l)==1:
            ind=l[0]
            return self.tps.element(ind)
        else:
            indtuple=l
            return self.tps.element(indtuple)


    def cst(self):
        return self.tps.cst()

    def linear(self):
        cdef PyDTPSA result
        result=PyDTPSA(0.0)
        result.tps=self.tps.linear()
        return result

    def __getitem__(self, unsigned long ind):
        return self.tps.element(ind)

    def __iadd__(PyDTPSA self, a):
        cdef PyDTPSA rhs
        if isinstance(a, PyDTPSA):
            rhs=<PyDTPSA>a
        elif isinstance(a, (float,int)):
            rhs=PyDTPSA(float(a))
        self.tps.add_to(rhs.tps)
        return self

    def __isub__(PyDTPSA self, a):
        cdef PyDTPSA rhs
        if isinstance(a, PyDTPSA):
            rhs=<PyDTPSA>a
        elif isinstance(a, (float,int)):
            rhs=PyDTPSA(float(a))
        self.tps.sub_to(rhs.tps)
        return self

    def __imul__(PyDTPSA self, a):
        cdef PyDTPSA rhs
        if isinstance(a, PyDTPSA):
            rhs=<PyDTPSA>a
        elif isinstance(a, (float,int)):
            rhs=PyDTPSA(float(a))
        self.tps.mul_to(rhs.tps)
        return self

    def __idiv__(PyDTPSA self, a):
        cdef PyDTPSA rhs
        if isinstance(a, PyDTPSA):
            rhs=<PyDTPSA>a
        elif isinstance(a, (float,int)):
            rhs=PyDTPSA(float(a))
        self.tps.div_to(rhs.tps)
        return self

    def __itruediv__(PyDTPSA self, a):
        cdef PyDTPSA rhs
        if isinstance(a, PyDTPSA):
            rhs=<PyDTPSA>a
        elif isinstance(a, (float,int)):
            rhs=PyDTPSA(float(a))
        self.tps.div_to(rhs.tps)
        return self


    def __add__(a,b):
        cdef PyDTPSA lhs
        if isinstance(a, PyDTPSA):
            lhs=PyDTPSA(copy_tpsa=<PyDTPSA> a)
            lhs.__iadd__(b)
        elif isinstance(b, PyDTPSA):
            print("here")
            lhs=PyDTPSA(copy_tpsa=<PyDTPSA> b)
            lhs.__iadd__(a)
        return lhs

    def __sub__(a,b):
        cdef PyDTPSA lhs
        if isinstance(a, PyDTPSA):
            lhs=PyDTPSA(copy_tpsa=<PyDTPSA> a)
            lhs.__isub__(b)
        elif isinstance(b, PyDTPSA):
            lhs=PyDTPSA(copy_tpsa=<PyDTPSA> b)
            lhs.__imul__(PyDTPSA(-1.0))
            lhs.__iadd__(a)
        return lhs

    def __mul__(a,b):
        cdef PyDTPSA lhs
        if isinstance(a, PyDTPSA):
            lhs=PyDTPSA(copy_tpsa=<PyDTPSA> a)
            lhs.__imul__(b)
        elif isinstance(b, PyDTPSA):
            lhs=PyDTPSA(copy_tpsa=<PyDTPSA> b)
            lhs.__imul__(a)
        return lhs

    def __truediv__(a,b):
        cdef PyDTPSA lhs
        if isinstance(a, PyDTPSA):
            lhs=PyDTPSA(copy_tpsa=<PyDTPSA> a)
            lhs.__itruediv__(b)
        elif isinstance(b, PyDTPSA):
            lhs=PyDTPSA(value=a)
            lhs.__itruediv__(b)
        return lhs

    def __div__(a,b):
        cdef PyDTPSA lhs
        if isinstance(a, PyDTPSA):
            lhs=PyDTPSA(copy_tpsa=<PyDTPSA> a)
            lhs.__idiv__(b)
        elif isinstance(b, PyDTPSA):
            lhs=PyDTPSA(a)
            lhs.__idiv__(b)
        return lhs

    def __repr__(self):
        return self.tps.print_to_string().decode("utf-8")


cdef class PyCTPSA:
    cdef tpsa_def.CTPS[cplx] tps

    def __cinit__(self, a=0.0, int dim_ind=0, PyCTPSA copy_tpsa=None, vector[cplx] input_map=[]):

        if copy_tpsa is not None:
            self.tps=copy_tpsa.tps
            return

        if len(input_map)>0:
            self.tps.set_map(input_map)
            return

        if tpsa_def.CTPS[cplx].Get_Max_Degree()<0:
            print("The class is not initialized yet.")

        if isinstance(a, (int, float, complex)):
            if dim_ind>0:
                if dim_ind <=tpsa_def.CTPS[cplx].Get_TPS_Dim():
                    self.tps.assign(complex(a), int(dim_ind))
                else:
                    raise IndexError ("TPSA dimension is out of range.")
            else:
                self.tps.assign(complex(a))



    @classmethod
    def initialize(cls, int dim, int max_order):
        tpsa_def.CTPS[cplx].Initialize(dim, max_order)

    @classmethod
    def get_power_index(cls, int ind):
        return tpsa_def.CTPS[double].Get_Power_Index(ind)

    @classmethod
    def get_max_degree(cls):
        return tpsa_def.CTPS[double].Get_Max_Degree()

    @classmethod
    def get_max_terms(cls):
        return tpsa_def.CTPS[double].Get_Max_Terms()


    def get_dim(self):
        return self.tps.get_dim()
    def get_degree(self):
        return self.tps.get_degree()
    def get_term(self):
        return self.tps.get_terms()


    def indices(self):
        cdef vector[cplx] tpsmap=self.tps.get_map()
        return tpsmap

    def find_index(self, vector[int] power_ind):
        cdef unsigned long temp
        temp=self.tps.find_index(power_ind)
        return temp

    def find_power(self, unsigned long n):
        cdef vector[int] temp = self.tps.find_power(n)
        return temp

    def evaluate(self, vector[cplx] values):
        return self.tps.evaluate(values)


    def composition(self, values):
        cdef vector[tpsa_def.CTPS[cplx]] vect
        cdef PyCTPSA a
        for a in values:
            vect.push_back(a.tps)
        cdef PyCTPSA result
        result=PyCTPSA(0.0)
        result.tps=self.tps.evaluate(vect)
        return result

    def derivative(self, int variable, int order):

        cdef PyCTPSA result
        result=PyCTPSA(0.0)
        result.tps=self.tps.derivative(variable, order)
        return result

    def integrate(self, int variable, cplx a0):
        cdef PyCTPSA result
        result=PyCTPSA(0.0)
        result.tps=self.tps.integrate(variable, a0)
        return result


    def conjugate(self, int mode=1):
        cdef PyCTPSA result
        result=PyCTPSA(0.0)
        result.tps=self.tps.conjugate(mode)
        return result

    def element(self, *l):
        cdef vector[int] indtuple
        cdef int ind
        if len(l)==1:
            ind=l[0]
            return self.tps.element(ind)
        else:
            indtuple=l
            return self.tps.element(indtuple)

    def cst(self):

        return self.tps.cst()

    def linear(self):
        cdef PyCTPSA result
        result=PyCTPSA(0.0)
        result.tps=self.tps.linear()
        return result

    def __getitem__(self, unsigned long ind):
        return self.tps.element(ind)

    def __iadd__(PyCTPSA self, a):
        cdef PyCTPSA rhs
        if isinstance(a, PyCTPSA):
            rhs=<PyCTPSA>a
        elif isinstance(a, (float,int, complex)):
            rhs=PyCTPSA(complex(a))
        self.tps.add_to(rhs.tps)
        return self

    def __isub__(PyCTPSA self, a):
        cdef PyCTPSA rhs
        if isinstance(a, PyCTPSA):
            rhs=<PyCTPSA>a
        elif isinstance(a, (float,int, complex)):
            rhs=PyCTPSA(complex(a))
        self.tps.sub_to(rhs.tps)
        return self

    def __imul__(PyCTPSA self, a):
        cdef PyCTPSA rhs
        if isinstance(a, PyCTPSA):
            rhs=<PyCTPSA>a
        elif isinstance(a, (float,int, complex)):
            rhs=PyCTPSA(complex(a))
        self.tps.mul_to(rhs.tps)
        return self

    def __idiv__(PyCTPSA self, a):
        cdef PyCTPSA rhs
        if isinstance(a, PyCTPSA):
            rhs=<PyCTPSA>a
        elif isinstance(a, (float,int, complex)):
            rhs=PyCTPSA(complex(a))
        self.tps.div_to(rhs.tps)
        return self

    def __itruediv__(PyCTPSA self, a):
        cdef PyCTPSA rhs
        if isinstance(a, PyCTPSA):
            rhs=<PyCTPSA>a
        elif isinstance(a, (float,int, complex)):
            rhs=PyCTPSA(complex(a))
        self.tps.div_to(rhs.tps)
        return self




    def __add__(a,b):
        cdef PyCTPSA lhs
        if isinstance(a, PyCTPSA):
            lhs=PyCTPSA(copy_tpsa=<PyCTPSA> a)
            lhs.__iadd__(b)
        elif isinstance(b, PyCTPSA):
            print("here")
            lhs=PyCTPSA(copy_tpsa=<PyCTPSA> b)
            lhs.__iadd__(a)
        return lhs

    def __sub__(a,b):
        cdef PyCTPSA lhs
        if isinstance(a, PyCTPSA):
            lhs=PyCTPSA(copy_tpsa=<PyCTPSA> a)
            lhs.__isub__(b)
        elif isinstance(b, PyCTPSA):
            lhs=PyCTPSA(copy_tpsa=<PyCTPSA> b)
            lhs.__imul__(PyCTPSA(-1.0))
            lhs.__iadd__(a)
        return lhs

    def __mul__(a,b):
        cdef PyCTPSA lhs
        if isinstance(a, PyCTPSA):
            lhs=PyCTPSA(copy_tpsa=<PyCTPSA> a)
            lhs.__imul__(b)
        elif isinstance(b, PyCTPSA):
            lhs=PyCTPSA(copy_tpsa=<PyCTPSA> b)
            lhs.__imul__(a)
        return lhs

    def __truediv__(a,b):
        cdef PyCTPSA lhs
        if isinstance(a, PyCTPSA):
            lhs=PyCTPSA(copy_tpsa=<PyCTPSA> a)
            lhs.__itruediv__(b)
        elif isinstance(b, PyCTPSA):
            lhs=PyCTPSA(a)
            lhs.__itruediv__(b)
        return lhs

    def __div__(a,b):
        cdef PyCTPSA lhs
        if isinstance(a, PyCTPSA):
            lhs=PyCTPSA(copy_tpsa=<PyCTPSA> a)
            lhs.__idiv__(b)
        elif isinstance(b, PyCTPSA):
            lhs=PyCTPSA(a)
            lhs.__idiv__(b)
        return lhs

    def __repr__(self):
        return self.tps.print_to_string().decode("utf-8")






def inv(a):
    cdef PyCTPSA resultc
    cdef PyDTPSA resultd
    if isinstance(a, PyDTPSA):
        resultd=PyDTPSA(0.0)
        resultd.tps=tpsa_def.inv((<PyDTPSA>a).tps)
        return resultd
    if isinstance(a, PyCTPSA):
        resultc=PyCTPSA(0.0)
        resultc.tps=tpsa_def.inv((<PyCTPSA>a).tps)
        return resultc

def exp(a):
    cdef PyCTPSA resultc
    cdef PyDTPSA resultd
    if isinstance(a, PyDTPSA):
        resultd=PyDTPSA(0.0)
        resultd.tps=tpsa_def.exp((<PyDTPSA>a).tps)
        return resultd
    if isinstance(a, PyCTPSA):
        resultc=PyCTPSA(0.0)
        resultc.tps=tpsa_def.exp((<PyCTPSA>a).tps)
        return resultc

def sqrt(a):
    cdef PyCTPSA resultc
    cdef PyDTPSA resultd
    if isinstance(a, PyDTPSA):
        if a.cst()<=0:
            raise ZeroDivisionError("The argument below zero in square root")
        resultd=PyDTPSA(0.0)
        resultd.tps=tpsa_def.sqrt((<PyDTPSA>a).tps)
        return resultd
    if isinstance(a, PyCTPSA):
        if abs(a.cst())==0:
            raise ZeroDivisionError("The argument zero in square root")
        resultc=PyCTPSA(0.0)
        resultc.tps=tpsa_def.sqrt((<PyCTPSA>a).tps)
        return resultc

def pow(a, float ord):
    cdef PyCTPSA resultc
    cdef PyDTPSA resultd
    if isinstance(a, PyDTPSA):
        resultd=PyDTPSA(0.0)
        resultd.tps=tpsa_def.pow((<PyDTPSA>a).tps, ord)
        return resultd
    if isinstance(a, PyCTPSA):
        resultc=PyCTPSA(0.0)
        resultc.tps=tpsa_def.pow((<PyCTPSA>a).tps, ord)
        return resultc

def log(a):
    cdef PyCTPSA resultc
    cdef PyDTPSA resultd
    if isinstance(a, PyDTPSA):
        if a.cst()<=0:
            raise ZeroDivisionError("The argument below zero in log")
        resultd=PyDTPSA(0.0)
        resultd.tps=tpsa_def.log((<PyDTPSA>a).tps)
        return resultd
    if isinstance(a, PyCTPSA):
        if abs(a.cst())==0:
            raise ZeroDivisionError("The argument zero in log")
        resultc=PyCTPSA(0.0)
        resultc.tps=tpsa_def.log((<PyCTPSA>a).tps)
        return resultc

def sin(a):
    cdef PyCTPSA resultc
    cdef PyDTPSA resultd
    if isinstance(a, PyDTPSA):
        resultd=PyDTPSA(0.0)
        resultd.tps=tpsa_def.sin((<PyDTPSA>a).tps)
        return resultd
    if isinstance(a, PyCTPSA):
        resultc=PyCTPSA(0.0)
        resultc.tps=tpsa_def.sin((<PyCTPSA>a).tps)
        return resultc

def arcsin(a):
    cdef PyCTPSA resultc
    cdef PyDTPSA resultd
    if isinstance(a, PyDTPSA):
        resultd=PyDTPSA(0.0)
        resultd.tps=tpsa_def.arcsin((<PyDTPSA>a).tps)
        return resultd
    if isinstance(a, PyCTPSA):
        resultc=PyCTPSA(0.0)
        resultc.tps=tpsa_def.arcsin((<PyCTPSA>a).tps)
        return resultc

def cos(a):
    cdef PyCTPSA resultc
    cdef PyDTPSA resultd
    if isinstance(a, PyDTPSA):
        resultd=PyDTPSA(0.0)
        resultd.tps=tpsa_def.cos((<PyDTPSA>a).tps)
        return resultd
    if isinstance(a, PyCTPSA):
        resultc=PyCTPSA(0.0)
        resultc.tps=tpsa_def.cos((<PyCTPSA>a).tps)
        return resultc

def arccos(a):
    cdef PyCTPSA resultc
    cdef PyDTPSA resultd
    if isinstance(a, PyDTPSA):
        resultd=PyDTPSA(0.0)
        resultd.tps=tpsa_def.arccos((<PyDTPSA>a).tps)
        return resultd
    if isinstance(a, PyCTPSA):
        resultc=PyCTPSA(0.0)
        resultc.tps=tpsa_def.arccos((<PyCTPSA>a).tps)
        return resultc

def tan(a):
    cdef PyCTPSA resultc
    cdef PyDTPSA resultd
    if isinstance(a, PyDTPSA):
        resultd=PyDTPSA(0.0)
        resultd.tps=tpsa_def.tan((<PyDTPSA>a).tps)
        return resultd
    if isinstance(a, PyCTPSA):
        resultc=PyCTPSA(0.0)
        resultc.tps=tpsa_def.tan((<PyCTPSA>a).tps)
        return resultc

def sinh(a):
    cdef PyCTPSA resultc
    cdef PyDTPSA resultd
    if isinstance(a, PyDTPSA):
        resultd=PyDTPSA(0.0)
        resultd.tps=tpsa_def.sinh((<PyDTPSA>a).tps)
        return resultd
    if isinstance(a, PyCTPSA):
        resultc=PyCTPSA(0.0)
        resultc.tps=tpsa_def.sinh((<PyCTPSA>a).tps)
        return resultc

def cosh(a):
    cdef PyCTPSA resultc
    cdef PyDTPSA resultd
    if isinstance(a, PyDTPSA):
        resultd=PyDTPSA(0.0)
        resultd.tps=tpsa_def.cosh((<PyDTPSA>a).tps)
        return resultd
    if isinstance(a, PyCTPSA):
        resultc=PyCTPSA(0.0)
        resultc.tps=tpsa_def.cosh((<PyCTPSA>a).tps)
        return resultc