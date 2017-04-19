from distutils.core import setup, Extension

tpsa_module = Extension('_tpsalib',
                       sources=['tpsalib_wrap.cxx', '../src/polymap.cpp', '../src/mathfunc.cpp'],
                       )

setup(name='tpsalib',
      version='0.1',
      author="Y. Hao",
      description="""Truncate Power Series Algorithm""",
      ext_modules=[tpsa_module],
      py_modules=["tpsalib"],
      )
