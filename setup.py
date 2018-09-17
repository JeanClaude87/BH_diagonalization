from distutils.core import setup
from Cython.Build import cythonize

setup(
	ext_modules = cythonize('bose.py', 'hamiltonian.py', 'function.py', 'observables.py')
)

# python setup.py build_ext --inplace
