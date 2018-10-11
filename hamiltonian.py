import numpy as np
from scipy.sparse import csc_matrix
from scipy.sparse import linalg as linalgS
from numpy import linalg as lin
import time
import function as ff
import hamiltonian_parity as ham_par
from joblib import Parallel, delayed

import multiprocessing

def bose_Hamiltonian (**args):

#Hamiltonian needs as input:

#...... Parameter: DIM_H

	DIM_H 	 = np.int(args.get("DIM_H"))
	ll 		 = args.get("ll")	
	nn 		 = args.get("nn")		

#Hamiltonian returns:
#...... Sparse or Dense matrix. By default is SPARSE
	
	mat_type = args.get("mat_type")
	if mat_type == None:
		mat_type = 'Sparse'


#...... Creates i,j,H[i,j]

	cores_num = np.int(args.get("cores_num"))

	if cores_num > 1:
		
		Hamiltonian = csc_matrix(([], ([], [])), shape=(DIM_H,DIM_H), dtype=np.double)
		
		step = DIM_H // cores_num
		split = Parallel(n_jobs=cores_num)(delayed(parallel_evaluate_ham)(step*k, step*(k+1),**args) for k in range(cores_num))

		for i in range(cores_num):
			
			Hamiltonian += split[i]

	else:

		X0 = [0]*DIM_H
		Y0 = [0]*DIM_H
		A0 = [0]*DIM_H

#............	
#............	JOAN PARALLELIZZA il LOOOOOOP a tope
#............

		pool = multiprocessing.Pool(2)
		X1,Y1,A1 = pool.map(functools.partial(evaluate_ham, **args), range(DIM_H))

		'''
		for i in range(DIM_H):

			X0[i],Y0[i],A0[i] = evaluate_ham(i, **args)
		'''

		#here we flatten the arrays
		X1 = [item for sublist in X0 for item in sublist]
		Y1 = [item for sublist in Y0 for item in sublist]
		A1 = [item for sublist in A0 for item in sublist]

		Hamiltonian = csc_matrix((A1, (X1,Y1)), shape=(DIM_H,DIM_H), dtype=np.double)

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

	if int_val != 0.0:
	
		#---- INTERACTION = we store i,i,int_val !!!!
		ham_ind1.append( i )
		ham_ind2.append( i )
		ham_val.append( int_val )

##----- KINETIC
	for hop in HOP_list:
		hop_state_bin = ff.TO_bin(state)^ff.TO_bin(hop)
		
		# we cut states with not N particles
		if ff.one_count(hop_state_bin) == nn:

			hop_state     = ff.TO_con(hop_state_bin,ll+nn-1)

			j = ff.get_index(hop_state,**args)	

			kin_val = t*action_hopping(i,j,**args)

			if kin_val != 0.0:

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
			
			if kin_val != 0.0:

	#---- KINETIC = we store i,j,kin_val !!!!
				ham_ind1.append( i )
				ham_ind2.append( j )
				ham_val.append( kin_val )

		if state[-1] == '1':

			PBC_newstate2 = state[-1]+state[:-1]
			j = ff.get_index(PBC_newstate2,**args)	

			kin_val = t*action_hopping(i,j,**args)
			
			if kin_val != 0.0:

	#---- KINETIC = we store i,j,kin_val !!!!
				ham_ind1.append( i )
				ham_ind2.append( j )
				ham_val.append( kin_val )

	return [ham_ind1, ham_ind2, ham_val]



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
	int_val = 0.5*U*np.dot(bosecon,bosecon-1.)

	return int_val

def diagonalization(Hamiltonian,**args):

	DIM_H    = args.get("DIM_H")
	mat_type = args.get("mat_type")
	parity	 = args.get("parity")
	num_eig	 = args.get("n_diag_state")

	if num_eig >= DIM_H:
		num_eig = DIM_H-2

	if mat_type == 'Sparse':

		A,B = linalgS.eigsh(Hamiltonian, k=num_eig, which='SA', return_eigenvectors=True)

	else:

		A,B   = lin.eigh(Hamiltonian)

	if parity == 'True':
		
		B = ham_par.vectors_parity_symmetrize(B,**args)

	return A,B





