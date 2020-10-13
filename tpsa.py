from .cython import tpsalib as tlib
import cmath

class tpsa(object):
    dimension = 0
    max_order = 0
    initialized = False
    string_repr=None
    def __init__(self, value=0.0, variable=0, dtype=float, tps=None):
        if tps is not None:
            if dtype in (float, complex):
                self.dtype=dtype
            else:
                raise ValueError("Unknown type")
            self._tps=tps;
            return
        if tpsa.initialized == False:
            print('TPSA class has to be initialized')
            exit(-1)

        if dtype == float:
            self._tps = tlib.PyDTPSA(float(value), variable)
            self.dtype=float

        elif dtype == complex:
            self._tps = tlib.PyCTPSA(complex(value), variable)
            self.dtype = complex

        else:
            raise ValueError("Unknown type")

    @classmethod
    def initialize(cls, dim, order):
        tpsa.dimension = dim
        tpsa.max_order = order
        tpsa.initialized = True
        tlib.PyDTPSA.initialize(dim, order)
        tlib.PyCTPSA.initialize(dim, order)

    @classmethod
    def set_variable_name(cls, string_rep):
        tpsa.string_repr = string_rep.split()
        if len(tpsa.string_repr) < tpsa.dimension:
            tpsa.string_repr = None
            print("Invalid number of representing setting")

    @classmethod
    def get_power_index(cls, ind):
        return tlib.PyDTPSA.get_power_index(ind)

    @classmethod
    def get_max_degree(cls):
        return tlib.PyDTPSA.get_max_degree()

    @classmethod
    def get_max_terms(cls):
        return tlib.PyDTPSA.get_max_terms()

    def cst(self):
        return self._tps.cst()

    '''@property
    def indices(self):
        power_list=[]
        value_list=[]
        for i in range(self.get_terms()):
            power = self.find_power(i)
            value = self.element(power)
            if value==0 and i>0:
                continue

            power_list.append(power[1:])
            value_list.append(value)
        return power_list, value_list
'''
    @property
    def indices(self):
        return self._tps.indices()

    def find_index(self, power_list):
        return self._tps.find_index(power_list)

    def find_power(self, index):
        return self._tps.find_power(index)

    def get_dim(self):
        return self._tps.get_dim()

    def get_degree(self):
        return self._tps.get_degree()

    def get_terms(self):
        return self._tps.get_term()

    def element(self, l):
        if isinstance(l,int):
            return self._tps.element(l)
        elif isinstance(l, list):
            return self._tps.element(*l)



    def derivative(self, dim, order=1):
        result=tpsa(0.0, dtype=self.dtype)
        result._tps=self._tps.derivative(dim, order)
        return result
    def integrate(self, dim, a0):
        result=tpsa(0.0, dtype=self.dtype)
        result._tps=self._tps.integrate(dim, a0)
        return result

    def __iadd__(self, other):
        if isinstance(other, tpsa):
            self._tps+=other._tps
        else:
            self = self + other
        return self

    def __isub__(self, other):
        if isinstance(other, tpsa):
            self._tps -= other._tps
        else:
            self = self - other
        return self

    def __imul__(self, other):
        if isinstance(other, tpsa):
            self._tps*=other._tps
        else:
            self = self * other
        return self

    def __idiv__(self, other):
        if isinstance(other, tpsa):
            self._tps /= other._tps
        else:
            self = self / other
        return self

    def __add__(self, other):
        result = tpsa(0.0, dtype=self.dtype)
        if isinstance(other, tpsa):
            result._tps = self._tps + other._tps
        else:
            result._tps = self._tps + other
        return result

    def __radd__(self, other):
        result = tpsa(0.0, dtype=self.dtype)
        result._tps = self._tps + other
        return result

    def __sub__(self, other):
        result = tpsa(0.0, dtype=self.dtype)
        if isinstance(other, tpsa):
            result._tps = self._tps - other._tps
        else:
            result._tps = self._tps - other
        return result

    def __rsub__(self, other):
        result = tpsa(0.0, dtype=self.dtype)
        result._tps = other - self._tps
        return result

    def __mul__(self, other):
        result = tpsa(0.0, dtype=self.dtype)
        if isinstance(other, tpsa):
            result._tps = self._tps * other._tps
        else:
            result._tps = self._tps * other
        return result

    def __rmul__(self, other):
        result = tpsa(0.0, dtype=self.dtype)
        result._tps = self._tps * other
        return result

    def __div__(self, other):
        result = tpsa(0.0, dtype=self.dtype)
        if isinstance(other, tpsa):
            result._tps = self._tps / other._tps
        else:
            result._tps = self._tps / other
        return result

    def __rdiv__(self, other):
        result = tpsa(0.0, dtype=self.dtype)
        result._tps = other / self._tps
        return result

    def __truediv__(self, other):
        result = tpsa(0.0, dtype=self.dtype)
        if isinstance(other, tpsa):
            result._tps = self._tps / other._tps
        else:
            result._tps = self._tps / other
        return result

    def __rtruediv__(self, other):
        result = tpsa(0.0, dtype=self.dtype)
        result._tps = other / self._tps
        return result

    def __neg__(self):
        return self * (-1.0)

    def __pos__(self):
        return self

    def __repr__(self):
        return self._tps.__repr__()

def inv(a):
    if isinstance(a, tpsa):
        return tpsa(tps=tlib.inv(a._tps), dtype=a.dtype)
    else:
        return 1.0/a

def exp(a):
    if isinstance(a, tpsa):
        return tpsa(tps=tlib.exp(a._tps), dtype=a.dtype)
    else:
        return cmath.exp(a)

def log(a):
    if isinstance(a, tpsa):
        return tpsa(tps=tlib.log(a._tps), dtype=a.dtype)
    else:
        return cmath.log(a)

def sqrt(a):
    if isinstance(a, tpsa):
        return tpsa(tps=tlib.sqrt(a._tps), dtype=a.dtype)
    else:
        return cmath.sqrt(a)

def pow(a, ind):
    if isinstance(a, tpsa):
        return tpsa(tps=tlib.pow(a._tps, ind), dtype=a.dtype)
    else:
        return None

def sin(a):
    if isinstance(a, tpsa):
        return tpsa(tps=tlib.sin(a._tps), dtype=a.dtype)
    else:
        return cmath.sin(a)
def arcsin(a):
    if isinstance(a, tpsa):
        return tpsa(tps=tlib.arcsin(a._tps), dtype=a.dtype)
    else:
        return cmath.asin(a)

def cos(a):
    if isinstance(a, tpsa):
        return tpsa(tps=tlib.cos(a._tps), dtype=a.dtype)
    else:
        return cmath.cos(a)
def arccos(a):
    if isinstance(a, tpsa):
        return tpsa(tps=tlib.arccos(a._tps), dtype=a.dtype)
    else:
        return cmath.acos(a)
def tan(a):
    if isinstance(a, tpsa):
        return tpsa(tps=tlib.tan(a._tps), dtype=a.dtype)
    else:
        return cmath.tan(a)

def sinh(a):
    if isinstance(a, tpsa):
        return tpsa(tps=tlib.sinh(a._tps), dtype=a.dtype)
    else:
        return cmath.sinh(a)
def cosh(a):
    if isinstance(a, tpsa):
        return tpsa(tps=tlib.cosh(a._tps), dtype=a.dtype)
    else:
        return cmath.cosh(a)

def initialize(dim, order, variable_name=None):
    tpsa.initialize(dim,order)
    if variable_name is not None:
        tpsa.set_variable_name(variable_name)

def get_dimension():
    return tpsa.dimension










