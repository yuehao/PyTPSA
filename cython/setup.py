from distutils.core import setup, Extension
from Cython.Build import cythonize

ext = Extension('tpsalib',
                sources=["tpsalib.pyx","../src/polymap.cpp", "../src/mathfunc.cpp"],  # additional source file(s)
                language="c++",             # generate C++ code
                #libraries=['gomp',],
                extra_compile_args=["-O3", "-std=c++11",
                                    #"-fopenmp",
                                    ],
                )

setup (name='tpsalib', ext_modules=cythonize(ext))

