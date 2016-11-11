from distutils.core import setup
from Cython.Build import cythonize

setup(
  name = '_fast_end_ramp_occ',
  ext_modules = cythonize("_fast_end_ramp_occ.pyx"),
)
