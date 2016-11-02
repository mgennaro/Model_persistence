from distutils.core import setup
from Cython.Build import cythonize

setup(
  name = 'cPixTrap',
  ext_modules = cythonize("cPixTraps.pyx")
    )
