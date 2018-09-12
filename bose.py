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


#........ LIST OF GLOBAL PARAMETERS

### -> ll, nn, tab_fact, DIM_H, BASE_bin, BASE_bose, CORR_BASE

ll           = 3
ff.ll        = ll

nn           = 2
ff.nn        = nn

ff.tab_fact  = tab_fact   = ff.fact_creation(nn+ll)

ff.DIM_H     = DIM_H      = ff.hilb_dim(nn,ll)

BASE_bin, BASE_bose       = ff.Base_prep()

ff.BASE_bin  = BASE_bin
ff.BASE_bose = BASE_bose

CORR_BASE    = ff.OUTER_creation(BASE_bose)
ff.CORR_BASE = CORR_BASE


#........ LIST OF HAMILTONIAN PARAMETERS

t=-1.
U=-2.0

BC=0




nstate = DIM_H

base_parity_ind = ff.base_parity()
print(base_parity_ind)

DIM_par_H = len(base_parity_ind)



ham_ind1, ham_ind2, ham_val = ff.bose_Hamiltonian_parity(base_parity_ind,BC,t,U,1)

print(ham_ind1)
print(ham_ind2)
#print(DIM_par_H)
#ED,VD = ff.diagonalization(ham_ind1, ham_ind2, ham_val, DIM_par_H, nstate, False)

#print(ED)

ham_ind1, ham_ind2, ham_val = ff.bose_Hamiltonian(BC,t,U)
ED,VD = ff.diagonalization(ham_ind1, ham_ind2, ham_val, DIM_H, nstate, False)

print(ED)


'''

ED,VD = ff.diagonalization(ham_ind1, ham_ind2, ham_val, nstate, False)
ES,VS = ff.diagonalization(ham_ind1, ham_ind2, ham_val, nstate, True)

print(ED)
print(ES)

print(VD)
print(VS)
'''

'''



st_ind = 5

for x in range(DIM_H-2):
	print(np.transpose(VD[x])
	print(VS[x])



for i in range(st_ind):

	print(E[i])

	corr0   = ff.NiNj   (   V[:,i], CORR_BASE)
	corr0_r = ff.NfixNr (5, V[:,i], CORR_BASE)

	dens0   = ff.density(   V[:,i])

	corr_name = str('NiNj')
	np.savetxt(PATH_now+os.sep+corr_name+str('_')+str(i)+str('.dat'), corr0, fmt='%.9f')

	corr_name_r = str('N5Nr')
	np.savetxt(PATH_now+os.sep+corr_name_r+str('_')+str(i)+str('.dat'), corr0_r, fmt='%.9f')

	dens_name = str('dens')
	np.savetxt(PATH_now+os.sep+dens_name+str('_')+str(i)+str('.dat'), dens0, fmt='%.9f')
'''
















