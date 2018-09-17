from distutils.core import setup
from Cython.Build import cythonize

setup(
	ext_modules = cythonize('function.py', 'observables.py', 'hamiltonian.py')
)

# python setup.py build_ext --inplace
