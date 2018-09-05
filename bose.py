import numpy as np
from math import factorial
import math
import itertools
from scipy.sparse import csc_matrix
from scipy.sparse import linalg

import function as ff
from numpy import linalg as LA

ll=14
nn=4

t=-1.
U=-1.

tab_fact = ff.fact_creation(nn+ll)

BC=0

BASE_num, BASE_bin = ff.Base_prep(ll,nn)
DIM_H = ff.hilb_dim(tab_fact,nn,ll)

ham_ind1,ham_ind2,ham_val = ff.bose_Hamiltonian(ll,nn,BC,t,U,BASE_bin,tab_fact)


numpy_ind1 = np.asarray(ham_ind1)
numpy_ind2 = np.asarray(ham_ind2)
numpy_val  = np.asarray(ham_val)


Hamiltonian = csc_matrix((numpy_val, (numpy_ind1, numpy_ind2)), shape=(DIM_H,DIM_H), dtype=np.double)
eig = linalg.eigsh(Hamiltonian, k=2, return_eigenvectors=False)
print(eig)

#uga = csc_matrix.todense(Hamiltonian)
#E,V = LA.eigh(uga)

#print(E)
















