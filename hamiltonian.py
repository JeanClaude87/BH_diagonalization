import numpy as np
from math import factorial
import math
import itertools
import os
from scipy.sparse import csc_matrix
from scipy.sparse import lil_matrix
from scipy.sparse import linalg as linalgS
from numpy import linalg as lin
from numpy import matlib
import time
import function as ff

from joblib import Parallel, delayed
import multiprocessing



def bose_Hamiltonian (**args):

#Hamiltonian needs as input:

#...... Parameter: DIM_H

	DIM_H 	 = args.get("DIM_H")

#Hamiltonian returns:
#...... Sparse or Dense matrix. By default is SPARSE
	
	mat_type = args.get("mat_type")
	if mat_type == None:
		mat_type = 'Sparse'


#...... Creates i,j,H[i,j]

	Hamiltonian = csc_matrix(([], ([], [])), shape=(DIM_H,DIM_H), dtype=np.double)
	
	num_cores = 2
	step = DIM_H // num_cores
	split = Parallel(n_jobs=num_cores)(delayed(parallel_evaluate_ham)(step*k, step*(k+1),**args) for k in range(num_cores))

	for i in range(len(split)):
		Hamiltonian += split[i]

	'''

	Hamiltonian = csc_matrix(([], ([], [])), shape=(DIM_H,DIM_H), dtype=np.double)

	for i in range(DIM_H):

		A, B, C = evaluate_ham(i,**args)

		Hamiltonian += csc_matrix((C, (A,B)), shape=(DIM_H,DIM_H), dtype=np.double)
	'''

	if mat_type == 'Dense':

		Hamiltonian = csc_matrix.todense(Hamiltonian)


	return Hamiltonian




def parallel_evaluate_ham (a,b,**args):
	
	DIM_H      = args.get("DIM_H")
	Hamiltonian = csc_matrix(([], ([], [])), shape=(DIM_H,DIM_H), dtype=np.double)
	
	for i in range(DIM_H):
		A, B, C = evaluate_ham(i,**args)
		Hamiltonian += csc_matrix((C, (A,B)), shape=(DIM_H,DIM_H), dtype=np.double)

	return Hamiltonian




def evaluate_ham(i,**args):

#...... Parameter: BC, t, U, DIM_H, nn, ll

	nn 		 = args.get("nn")
	ll 		 = args.get("ll")	
	BC 		 = args.get("BC")
	t 		 = args.get("t")	
	U 		 = args.get("U")
	DIM_H 	 = args.get("DIM_H")	
	tab_fact = args.get("tab_fact")

#...... Functions in ham.py: action_hopping, action_interactions
#...... Fundamental tables:  BASE_bin, HOP_list

	BASE_bin = args.get("BASE_bin")
	HOP_list = args.get("HOP_list")	

#Hamiltonian returns:
#...... Sparse or Dense matrix. By default is SPARSE
	
	mat_type = args.get("mat_type")
	if mat_type == None:
		mat_type = 'Sparse'

	ham_ind1 = []
	ham_ind2 = []
	ham_val  = []

	state = BASE_bin[i]

##----- INTERACTIONS

	int_val = action_interactions(state,**args)

	#---- INTERACTION = we store i,i,int_val !!!!
	ham_ind1.append( i )
	ham_ind2.append( i )
	ham_val.append( int_val )

##----- KINETIC
	for hop in HOP_list:
		hop_state_bin = ff.TO_bin(state)^ff.TO_bin(hop)
		hop_state     = ff.TO_con(hop_state_bin,ll+nn-1)
		
		# we cut states with not N particles
		if ff.one_count(hop_state_bin) == nn:

			j = ff.get_index(hop_state,**args)	
			kin_val = t*action_hopping(i,j,**args)

	#---- KINETIC = we store i,j,kin_val !!!!
			ham_ind1.append( i )
			ham_ind2.append( j )
			ham_val.append( kin_val )

	# here we put the PERIODIC
	#...................BC=0 periodic, BC=1 open

	if BC == 0:
		if state[0] == '1':

			PBC_newstate1 = state[1:]+state[0]				
			j = ff.get_index(PBC_newstate1,**args)	

			kin_val = t*action_hopping(i,j,**args)
			
	#---- KINETIC = we store i,j,kin_val !!!!
			ham_ind1.append( i )
			ham_ind2.append( j )
			ham_val.append( kin_val )

		if state[-1] == '1':

			PBC_newstate2 = state[-1]+state[:-1]
			j = ff.get_index(PBC_newstate2,**args)	

			kin_val = t*action_hopping(i,j,**args)
			
	#---- KINETIC = we store i,j,kin_val !!!!
			ham_ind1.append( i )
			ham_ind2.append( j )
			ham_val.append( kin_val )

	return ham_ind1, ham_ind2, ham_val



def action_hopping(x,y,**args):

	ll 		 = args.get("ll")	
	BASE_bin = args.get("BASE_bin")


	state_x   = BASE_bin[x]
	bosecon_x = ff.TO_bose_conf(state_x,ll)

	state_y   = BASE_bin[y]
	bosecon_y = ff.TO_bose_conf(state_y,ll)

	jump_c  =	np.argmax(bosecon_x-bosecon_y)
	jump_cd  =	np.argmin(bosecon_x-bosecon_y)

	result	  = np.sqrt((bosecon_x[jump_cd]+1)*(bosecon_x[jump_c])) 
	
	return result


def action_interactions(state,**args):


	ll 		 = args.get("ll")
	U 		 = args.get("U")

	bosecon = ff.TO_bose_conf(state,ll)

	int_val = 0
	for x in range(ll):
		nx = bosecon[x]
		int_val += U*nx*(nx-1.)/(2.) #+10**-3*np.random.random()

	return int_val

def diagonalization(Hamiltonian,num_eig,**args):

# num_eig 	-> how many eigenvalues: less than DIM_H


	DIM_H    = args.get("DIM_H")
	mat_type = args.get("mat_type")

	if num_eig >= DIM_H:
		num_eig = DIM_H-2

	if mat_type == 'Sparse':

		eig = linalgS.eigsh(Hamiltonian, k=num_eig, which='SA', return_eigenvectors=True)

	else:

		hamdens = csc_matrix.todense(Hamiltonian)
		eig  = lin.eigh(hamdens)
		eigl = list(eig)
		eigl[1] = np.squeeze(np.asarray(eig[1]))
		eig = eigl

	return eig
