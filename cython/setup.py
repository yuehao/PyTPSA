from distutils.core import setup, Extension
from Cython.Build import cythonize

ext = Extension('tpsalib',
                sources=["tpsalib.pyx","../src/polymap.cpp", "../src/mathfunc.cpp"],  # additional source file(s)
                language="c++",             # generate C++ code
                extra_compile_args=["-O3"],
                )

setup (name='tpsalib', ext_modules=cythonize(ext))

