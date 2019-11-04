import numpy as np
import scipy as sp
from scipy.sparse import linalg as linalgS
from scipy import linalg as lin
import function as ff
import hamiltonian_parity as ham_par


def evaluate_ham(i,**args):

#...... Parameter: BC, t, U, DIM_H, nn, ll

	nn 		 = args.get("nn")
	ll 		 = args.get("ll")	
	BC 		 = args.get("BC")

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

##----- Barrier

	int_bar = action_potential(state,**args)

	if int_bar != 0.0:
	
		#---- INTERACTION = we store i,i,int_val !!!!
		ham_ind1.append( i )
		ham_ind2.append( i )
		ham_val.append( int_bar )

##----- KINETIC_first
	for hop in HOP_list:
		hop_state_bin = ff.TO_bin(state)^ff.TO_bin(hop)
		
		# we cut states with not N particles
		if ff.one_count(hop_state_bin) == nn:

			hop_state     = ff.TO_con(hop_state_bin,ll+nn-1)

			j = ff.get_index(hop_state,**args)	
			kin_val = action_hopping(i,j,**args)

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
			kin_val = action_hopping(i,j,**args)
			
			if kin_val != 0.0:

	#---- KINETIC = we store i,j,kin_val !!!!
				ham_ind1.append( i )
				ham_ind2.append( j )
				ham_val.append( kin_val )

		if state[-1] == '1':

			PBC_newstate2 = state[-1]+state[:-1]
			j = ff.get_index(PBC_newstate2,**args)	
			kin_val = action_hopping(i,j,**args)
			
			if kin_val != 0.0:

	#---- KINETIC = we store i,j,kin_val !!!!
				ham_ind1.append( i )
				ham_ind2.append( j )
				ham_val.append( kin_val )

	return [ham_ind1, ham_ind2, ham_val]

def action_hopping(x,y,**args):

	ll 		 = args.get("ll")	
	BASE_bin = args.get("BASE_bin")
	t 		 = args.get("t")

	state_x   = BASE_bin[x]
	bosecon_x = ff.TO_bose_conf(state_x,ll)

	state_y   = BASE_bin[y]
	bosecon_y = ff.TO_bose_conf(state_y,ll)


	jump_cd  =	np.argmin(bosecon_x-bosecon_y)
	jump_c  =	np.argmax(bosecon_x-bosecon_y)

	rem = (jump_c+1)%ll

	if   rem == jump_cd:
		direction = 'right'
		hop       = np.conjugate(t) 
	elif jump_c == jump_cd:
		direction = 'stand'
	else:
		direction = 'left'
		hop       = t		

#	print(bosecon_x,bosecon_y,direction)

	result	  = np.sqrt((bosecon_x[jump_cd]+1)*(bosecon_x[jump_c])) 

	return hop*result


def action_hopping_second(x,y,**args):

	ll 		 = args.get("ll")	
	BASE_bin = args.get("BASE_bin")
	t 		 = args.get("t")

	state_x   = BASE_bin[x]
	bosecon_x = ff.TO_bose_conf(state_x,ll)

	state_y   = BASE_bin[y]
	bosecon_y = ff.TO_bose_conf(state_y,ll)


	jump_cd  =	np.argmin(bosecon_x-bosecon_y)
	jump_c  =	np.argmax(bosecon_x-bosecon_y)

	rem = (jump_c+1)%ll

	if   rem == jump_cd:
		direction = 'right'
		hop       = np.conjugate(t) 
	elif jump_c == jump_cd:
		direction = 'stand'
	else:
		direction = 'left'
		hop       = t		

#	print(bosecon_x,bosecon_y,direction)

	result	  = np.sqrt((bosecon_x[jump_cd]+1)*(bosecon_x[jump_c])) 

	return hop*result


def action_interactions(state,**args):

	ll 		 = args.get("ll")
	U 		 = args.get("U")

	bosecon = ff.TO_bose_conf(state,ll)
	int_val = 0.5*U*np.dot(bosecon,bosecon-1.)

	return int_val

def action_potential(state,**args):

	ll 		= args.get("ll")
	bar		= args.get("bar")	

	pot 	= np.zeros(ll, dtype=np.float)
	pot[0]  = bar

	bosecon = ff.TO_bose_conf(state,ll)
	bar_val = np.dot(pot,bosecon)
#	print(pot)
#	print(bosecon)
#	print(bar_val)

	return bar_val

def diagonalization(Hamiltonian,**args):

	DIM_H    = args.get("DIM_H")
	mat_type = args.get("mat_type")
	parity	 = args.get("parity")
	num_eig	 = args.get("n_diag_state")

	if num_eig <= 0:
		num_eig = 1
	
	if isinstance( Hamiltonian, sp.sparse.csc.csc_matrix):	
		if num_eig >= DIM_H-1:
			
			k0 = int(DIM_H/2)
			k1 = DIM_H - k0

			A0,B0   =  linalgS.eigsh( Hamiltonian, k=k0, which='SA', return_eigenvectors=True)
			A1m,B1m =  linalgS.eigsh(-Hamiltonian, k=k1, which='SA', return_eigenvectors=True)

			A1=-A1m
			B1=-B1m

			A = np.concatenate((A0,A1))
			B = np.concatenate((B0.T,B1.T)).T

		else:
			
			A,B = linalgS.eigsh(Hamiltonian, k=num_eig, which='SA', return_eigenvectors=True)

	else:

		if num_eig >= DIM_H:
			num_eig = DIM_H

		A,B = sp.linalg.eigh(Hamiltonian)

	if parity == 'True':
		
		B = ham_par.vectors_parity_symmetrize(B,**args)

	B = B[:, A.argsort()]
	A = np.sort(A)

	return A[:num_eig],B[:,:num_eig]





