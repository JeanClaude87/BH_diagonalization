import os

import numpy as np
from math import factorial
import math
import itertools
from scipy.sparse import csc_matrix
from scipy.sparse import linalg
from numpy import linalg as LA

import function as ff




corr_1b = False
dens 	= False



PATH_now = os.path.abspath('.')

ll=3
nn=2

tab_fact = ff.fact_creation(nn+ll)

t=-1.
U=-2.0

BC=0


BASE_bin, BASE_bose = ff.Base_prep(ll,nn)
DIM_H = ff.hilb_dim(tab_fact,nn,ll)

nstate = DIM_H

ham_ind1, ham_ind2, ham_val = ff.bose_Hamiltonian(ll,nn,BC,t,U,BASE_bin,tab_fact)

ED,VD = ff.diagonalization(ham_ind1, ham_ind2, ham_val, DIM_H, nstate, False)
ES,VS = ff.diagonalization(ham_ind1, ham_ind2, ham_val, DIM_H, nstate, True)

print(ED)
print(ES)

print(VD)
print(VS)


'''

CORR_BASE = ff.OUTER_creation(BASE_bose)

st_ind = 5

for x in range(DIM_H-2):
	print(np.transpose(VD[x])
	print(VS[x])



for i in range(st_ind):

	print(E[i])

	corr0   = ff.NiNj   (   V[:,i], CORR_BASE)
	corr0_r = ff.NfixNr (5, V[:,i], CORR_BASE)

	dens0   = ff.density(   V[:,i], BASE_bose)

	corr_name = str('NiNj')
	np.savetxt(PATH_now+os.sep+corr_name+str('_')+str(i)+str('.dat'), corr0, fmt='%.9f')

	corr_name_r = str('N5Nr')
	np.savetxt(PATH_now+os.sep+corr_name_r+str('_')+str(i)+str('.dat'), corr0_r, fmt='%.9f')

	dens_name = str('dens')
	np.savetxt(PATH_now+os.sep+dens_name+str('_')+str(i)+str('.dat'), dens0, fmt='%.9f')
'''
















