from distutils.core import setup
from Cython.Build import cythonize

setup(
#	ext_modules = cythonize( ['hamiltonian.py', 'hamiltonian_parity.py', 'function.py', 'observables.py', 'time_evolution.py'])
ext_modules = cythonize( ['observables.py','function.py' ])

)

# python3 setup.py build_ext --inplace
