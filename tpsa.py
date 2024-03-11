from .cython import tpsalib as tlib
import cmath




class tpsa(object):
    dimension = 0
    max_order = 0
    initialized = False
    def_string_repr=None
    def __init__(self, value=0.0, variable=0, dtype=float, tps=None, input_map=[]):
        self.string_repr=tpsa.def_string_repr
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

        if len(input_map) > 0:
            if dtype==float:
                self._tps=tlib.PyDTPSA(input_map=input_map)
                self.dtype=dtype
                return

            elif dtype==complex:
                self._tps=tlib.PyCTPSA(input_map=input_map)
                self.dtype = dtype
                return
            else:
                raise ValueError("Unknown type")

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
    def linear(self):
        result = tpsa(0.0, dtype=self.dtype)
        result._tps=self._tps.linear()
        return result

    def pvl(self):
        power_list=[]
        value_list=[]
        maxind=len(self.indices)
        for i in range(self.get_max_terms()):
            power = self.find_power(i)
            if i<maxind:
                value = self.element(power)
            else:
                value = 0
            #if value==0 and i>0:
            #    continue

            power_list.append(power[1:])
            value_list.append(value)
        return power_list, value_list


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

    def evaluate(self, values):
        if len(values)==tpsa.dimension and isinstance(values[0], self.dtype):
            return self._tps.evaluate(values)
        else:
            print("wrong vector length or wrong datatype")

    def composite(self, values):
        if len(values)==tpsa.dimension:
            try:
                tps_list=[v._tps for v in values]
            except:
                tps_list=[]
                for v in values:
                    if isinstance(v, tpsa):
                        tps_list.append(v._tps)
                    else:
                        tps_list.append(tpsa(v,dtype=self.dtype)._tps)

            temp = self._tps.composition(tps_list)
            return tpsa(tps=temp, dtype=self.dtype)
        else:
            print("wrong vector length, should be {}".format(tpsa.dimension))

    def derivative(self, dim, order=1):
        result=tpsa(0.0, dtype=self.dtype)
        result._tps=self._tps.derivative(dim, order)
        return result
    def integrate(self, dim, a0):
        result=tpsa(0.0, dtype=self.dtype)
        result._tps=self._tps.integrate(dim, a0)
        return result

    def conjugate(self, mode=''):
        # mode=1; z1, z1*, z2, z2* ...
        # mode=2; z1, z2,...,z1*, z2* ...
        modenum = 0
        if mode.upper() == 'COMPLEXPAIR' or mode.upper() == 'CP':
            modenum = 1
        elif mode.upper() == 'REAL' or mode.upper() == 'R':
            modenum = 2
        else:
            print("INVALID MODE, supported mode are 'ComplexPair' or 'Real'")
            return None
        result = tpsa(0.0, dtype=self.dtype)
        result._tps = self._tps.conjugate(modenum)
        return result

    def copy(self):
        return tpsa(dtype=self.dtype, tps=self._tps)


    def __iadd__(self, other):
        if isinstance(other, tpsa):
            self._tps += other._tps
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
        if isinstance(other, tpsa):
            result._tps = self._tps / other._tps
        else:
            result._tps = (tpsa(other, dtype=self.dtype))._tps / self._tps
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

def save(filename, *args):

    import numpy as np
    if len(args)==0:
        print('Nothing saved.')
        return
    elif len(args)==1 and isinstance(args[0],list):
        save(filename, *(args[0]))
        return

    with open(filename, 'wb') as f:
        np.save(f, np.array([len(args), args[0].dimension, args[0].max_order]))
        for arg in args:
            np.save(f, np.array(arg.indices))



def load(filename):
    import numpy as np
    ret=[]
    maxord=0
    maxterm=-1
    with open(filename, 'rb') as f:
        temp=np.load(f)
        print("Loading {} TPS with {} variables upto {} orders".format(temp[0],temp[1],temp[2]))
        if tpsa.initialized:
            if tpsa.dimension!=temp[1]:
                print('The saved TPSA had different dimension than current setting. Abort')
                return
            if tpsa.max_order<temp[2]:
                from scipy.special import comb
                maxterm=comb(tpsa.max_order+tpsa.dimension, tpsa.dimension, exact=True)
        else:
            tpsa.initialize(temp[1],temp[2])
        for i in range(temp[0]):
            temptps=np.load(f)
            if maxterm>0 and len(temp)>maxterm:
                temptps=temptps[0:maxterm]
                print('Warning, the tpsa is truncated from order {} to {}'.format(temp[2], tpsa.max_order))
            if isinstance(temptps[0], np.floating):
                dtype = float
            elif isinstance(temptps[0], np.complexfloating):
                dtype = complex
            else:
                raise ValueError("Unknown type")
            ret.append(tpsa(input_map=temptps.tolist(), dtype=dtype) )
        print("Success!")
        return ret


def inverse_map(list_of_tps, initial_trial=None, iteration_limit=0):
    import numpy as np
    if tpsa.dimension!=len(list_of_tps):
        print("The input is either over- or under- determined, inversion is not possible")
        return
    linear_map=np.eye(tpsa.dimension, dtype=list_of_tps[0].dtype)
    nlm=[]
    for i in range(tpsa.dimension):
        linear_map[i, :] = list_of_tps[i].linear().indices[1:]
        nlm.append(list_of_tps[i]-list_of_tps[i].linear())
    try:
        lminv=np.linalg.inv(linear_map)
    except:
        print("Linear map is non inversible, Abort")
        return None

    Iv = [tpsa(0, i + 1, dtype=list_of_tps[0].dtype) for i in range(tpsa.dimension)]
    results = [tpsa(0, i + 1, dtype=list_of_tps[0].dtype) for i in range(tpsa.dimension)]
    if initial_trial is not None:
        results=initial_trial
    temp = [tpsa(0, i + 1, dtype=list_of_tps[0].dtype) for i in range(tpsa.dimension)]
    iterations=tpsa.max_order
    if iteration_limit>0 and iteration_limit<iterations:
        iterations=iteration_limit
    for k in range(iterations):
        for i in range(tpsa.dimension):
            temp[i] = (Iv[i] - nlm[i].composite(results))
        for i in range(tpsa.dimension):
            results=lminv.dot(temp).tolist()
    return results