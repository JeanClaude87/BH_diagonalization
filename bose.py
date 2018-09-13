import os

import numpy as np
from math import factorial
import math
import itertools
from scipy.sparse import csc_matrix
from scipy.sparse import linalg
from numpy import linalg as LA
import time
import function as ff

np.set_printoptions(precision=4)

start = time.time()

#........ LIST OF GLOBAL PARAMETERS

### -> ll, nn, tab_fact, DIM_H, BASE_bin, BASE_bose, CORR_BASE, PATH_now

ll           = 30
ff.ll        = ll

nn           = 5
ff.nn        = nn

ff.tab_fact  = tab_fact   = ff.fact_creation(nn+ll)

ff.DIM_H     = DIM_H      = ff.hilb_dim(nn,ll)

print('dim H', DIM_H)

BASE_bin, BASE_bose       = ff.Base_prep()

ff.BASE_bin  = BASE_bin
ff.BASE_bose = BASE_bose

CORR_BASE    = ff.OUTER_creation(BASE_bose)
ff.CORR_BASE = CORR_BASE

PATH_now = os.path.abspath('.')
#ff.PATH_now


#........ LIST OF HAMILTONIAN PARAMETERS

t=-1.
U=-1.0

BC=0



nstate = DIM_H

base_parity_ind = ff.base_parity()

DIM_par_H = len(base_parity_ind)

end = time.time()
tempotras = (end - start)
print('preparation done', tempotras)

start = time.time()

ham_ind1, ham_ind2, ham_val = ff.bose_Hamiltonian(BC,t,U)

end = time.time()
tempotras = (end - start)
print ('Hamiltonian built', tempotras)

start = time.time()

Sp_Hamiltonian = ff.make_sparse_mat(ham_ind1, ham_ind2, ham_val, DIM_H)

#E,V   = ff.bose_Hamiltonian_parity(Sp_Hamiltonian,base_parity_ind)

end = time.time()
tempotras = (end - start)
print ('Hamiltonian sparse', tempotras)

start = time.time()

ED,VD = ff.diagonalization(Sp_Hamiltonian, 2, True)

end = time.time()
tempotras = (end - start)
print ('Hamiltonian diagonal', tempotras)

number_state = 2

for i in range(number_state):
	
	dens   = ff.density(VD[:,i])
	print(i, 'dens', dens)



'''
#	corr0   = ff.NiNj   (V[:,i])
#	corr0_r = ff.NfixNr (V[:,i],int(xx))

#	corr_name = str('NiNj')
#	np.savetxt(PATH_now+os.sep+corr_name+str('_')+str(i)+str('.dat'),   corr0,   fmt='%.3f')

#	corr_name_r = str('N_'+str(xx)+'Nr')
#	np.savetxt(PATH_now+os.sep+corr_name_r+str('_')+str(i)+str('.dat'), corr0_r, fmt='%.3f')

#	dens_name = str('dens')
#	np.savetxt(PATH_now+os.sep+dens_name+str('_')+str(i)+str('.dat'),   dens0,   fmt='%.3f')






#ham_ind1, ham_ind2, ham_val = ff.bose_Hamiltonian_parity(base_parity_ind,BC,t,U,1)

#print(ham_ind1)
#print(ham_ind2)
#print(DIM_par_H)
#ED,VD = ff.diagonalization(ham_ind1, ham_ind2, ham_val, DIM_par_H, nstate, False)

#print(ED)

#ham_ind1, ham_ind2, ham_val = ff.bose_Hamiltonian(BC,t,U)
#Sp_Hamiltonian = ff.make_sparse_mat(ham_ind1, ham_ind2, ham_val, DIM_H)

#ED,VD = ff.diagonalization(Sp_Hamiltonian, nstate, False)

#ff.print_hamiltonian(Sp_Hamiltonian)




ED,VD = ff.diagonalization(ham_ind1, ham_ind2, ham_val, nstate, False)
ES,VS = ff.diagonalization(ham_ind1, ham_ind2, ham_val, nstate, True)

print(ED)
print(ES)

print(VD)
print(VS)




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
















