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
import function_jp as ffjp

np.set_printoptions(precision=2)
PATH_now = os.path.abspath('.')


#........ LIST OF GLOBAL PARAMETERS

### -> ll, nn, tab_fact, DIM_H, BASE_bin, BASE_bose, CORR_BASE

ll           = 8
ff.ll        = ll

nn           = 6
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
U=-3.0

BC=0




nstate = DIM_H

base_parity_ind = ff.base_parity()

#~ print(base_parity_ind)

DIM_par_H = len(base_parity_ind)



################################################ Original Hamiltonian

print('Preparing Base')
time0=time.clock()
ham_ind1, ham_ind2, ham_val = ff.bose_Hamiltonian(BC,t,U)
time1=time.clock()
print ('')
print ('Time=',time0,time1,"%e" % (time1-time0))
print ('')


#~ print('original bose base',BASE_bose)
Sp_Hamiltonian = ff.make_sparse_mat(ham_ind1, ham_ind2, ham_val, DIM_H)

#~ ff.print_hamiltonian(Sp_Hamiltonian)

print('Diagonalizing')
time0=time.clock()
ED1 ,EV1 = ff.diagonalization(Sp_Hamiltonian,DIM_H-1,True)
time1=time.clock()
print ('')
print ('Time=',time0,time1,"%e" % (time1-time0))
print ('')

#~ print(ED1)


################################################ PIERO STUFF

print('Preparing Base')
time0=time.clock()
H_par = ff.bose_Hamiltonian_parity(Sp_Hamiltonian,base_parity_ind,BC,t,U,1)
time1=time.clock()
print ('')
print ('Time=',time0,time1,"%e" % (time1-time0))
print ('')

#~ ff.print_hamiltonian(H_par)

print('Diagonalizing')
time0=time.clock()
ED ,EV = ff.diagonalization(H_par,DIM_H-1,True)
time1=time.clock()
print ('')
print ('Time=',time0,time1,"%e" % (time1-time0))
print ('')


#~ print(ED)



################################################ My stuff

print('Preparing Base')
time0=time.clock()
ham_p_ind1, ham_p_ind2, ham_p_val = ffjp.bose_Hamiltonian_parity2(base_parity_ind,ham_ind1,ham_ind2,ham_val,BC,t,U,1)
time1=time.clock()
print ('')
print ('Time=',time0,time1,"%e" % (time1-time0))
print ('')


SP_Hamiltonian = ff.make_sparse_mat(ham_p_ind1, ham_p_ind2, ham_p_val, DIM_H)

#~ ff.print_hamiltonian(SP_Hamiltonian)

print('Diagonalizing')
time0=time.clock()
ED ,EV = ff.diagonalization(SP_Hamiltonian,2,False)
time1=time.clock()
print ('')
print ('Time=',time0,time1,"%e" % (time1-time0))
print ('')

#~ print(ED)









#~ difference = []
#~ difference_temp = ED[:-1]-ED1
#~ for i in difference_temp:
	#~ if i < 1e-10:
		#~ difference.append(0)
	#~ else:
		#~ difference.append(i)
#~ print ('difference with respect to sparse', difference)





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
















