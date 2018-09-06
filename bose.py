import os

import numpy as np
from math import factorial
import math
import itertools
from scipy.sparse import csc_matrix
from scipy.sparse import linalg
from numpy import linalg as LA

import function as ff


PATH_now = os.path.abspath('.')

ll=10
nn=2

tab_fact = ff.fact_creation(nn+ll)

t=-1.
U=-2.

BC=0

nstate = 5



BASE_bin, BASE_bose = ff.Base_prep(ll,nn)

DIM_H = ff.hilb_dim(tab_fact,nn,ll)

ham_ind1, ham_ind2, ham_val = ff.bose_Hamiltonian(ll,nn,BC,t,U,BASE_bin,tab_fact)

E,V = ff.diagonalization(ham_ind1,ham_ind2,ham_val,DIM_H, nstate)



CORR_BASE = ff.OUTER_creation(BASE_bose)

corr0 = ff.NiNj(V[:,0],CORR_BASE)
dens0 = ff.density(V[:,0],BASE_bose)

ext = str('.dat')

corr_name = str('NiNj')
np.savetxt(PATH_now+os.sep+corr_name+ext, corr0, fmt='%.9f')

dens_name = str('dens')
np.savetxt(PATH_now+os.sep+dens_name+ext, dens0, fmt='%.9f')

















