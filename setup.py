from distutils.core import setup
from Cython.Build import cythonize

setup(ext_modules=cythonize(['spindex4.pyx']))

# python setup.py build_ext --inplace