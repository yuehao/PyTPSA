from .lib import tpsalib as tlib
import copy

class tpsa(object):
    dimension = 0
    max_order = 0
    initialized = False
    def __init__(self, value=0.0, variable=0, dtype=float, copy_tps=None):
        if tpsa.initialized == False:
            print('TPSA class has to be initialized')
            exit(-1)

        if copy_tps is not None:
            self.__tps=copy_tps
            return;

        if dtype == float:
            self.__tps = tlib.DCTPS(value)
            if variable > 0 and variable <= tpsa.dimension:
                self.__tps.assign(value, variable)

        elif dtype == complex:
            self.__tps = tlib.CCTPS(value)
            if variable > 0 and variable <= tpsa.dimension:
                self.__tps.assign(value, variable)

        else:
            print("Unknown type")
            exit(-1)

    @classmethod
    def initialize(cls, dim, order):
        tpsa.dimension = dim
        tpsa.max_order = order
        tpsa.initialized = True
        tlib.DCTPS.Initialize(dim, order)
        tlib.CCTPS.Initialize(dim, order)

    def cst(self):
        return self.__tps.cst()

    def findindex(self, power_list):
        return self.__tps.findindex(power_list)

    def findpower(self, index):
        return self.__tps.findpower(index)

    def get_dim(self):
        return self.__tps.get_dim()

    def get_degree(self):
        return self.__tps.get_degree()

    def get_terms(self):
        return self.__tps.get_terms()

    def element(self, index):
        return self.__tps.element(index)

    def evaluate(self, vectors):
        return self.__tps.evaluate(vectors)

    def derivative(self, dim, order=1):
        return tpsa(copy_tps=self.__tps.derivative(dim, order))

    def __add__(self, other):
        return tpsa(copy_tps=(self.__tps + other))

    def __radd__(self, other):
        return tpsa(copy_tps=(self.__tps + other))

    def __sub__(self, other):
        return tpsa(copy_tps=(self.__tps - other))

    def __rsub__(self, other):
        return tpsa(copy_tps=(other - self.__tps))

    def __mul__(self, other):
        return tpsa(copy_tps=(self.__tps * other))

    def __rmul__(self, other):
        return tpsa(copy_tps=(self.__tps * other))

    def __div__(self, other):
        return tpsa(copy_tps=(self.__tps / other))

    def __rdiv__(self, other):
        return tpsa(copy_tps=(other / self.__tps))

    def __neg__(self):
        return tpsa(copy_tps=self.__tps.__neg__())

    def __pos__(self):
        return tpsa(copy_tps=self.__tps.__pos__())

    def __repr__(self):
        return self.__tps.__repr__()

def inv(a):
    return tpsa(copy_tps=tlib.inv(a.__tps))

def exp(a):
    return tpsa(copy_tps=tlib.exp(a.__tps))

def log(a):
    return tpsa(copy_tps=tlib.log(a.__tps))

def sqrt(a):
    return tpsa(copy_tps=tlib.sqrt(a.__tps))

def pow(a):
    return tpsa(copy_tps=tlib.pow(a.__tps))

def sin(a):
    return tpsa(copy_tps=tlib.sin(a.__tps))

def cos(a):
    return tpsa(copy_tps=tlib.cos(a.__tps))

def tan(a):
    return tpsa(copy_tps=tlib.tan(a.__tps))

def sinh(a):
    return tpsa(copy_tps=tlib.sinh(a.__tps))

def cosh(a):
    return tpsa(copy_tps=tlib.cosh(a.__tps))









