import os

import numpy as np
from numpy import matlib
from math import factorial
import math
import itertools
from scipy.sparse import csc_matrix
from scipy.sparse import linalg
from numpy import linalg as LA
import time

import hamiltonian        as ham
import function           as ff
import observables        as ob

#..................................parity transformation
# A states
def parity(state,**args):
	
	parity_state 	= state[::-1]

	index_ps = ff.get_index(parity_state,**args)
	
	return parity_state,index_ps


def base_parity(**args):

	DIM_H 	 = np.int(args.get("DIM_H"))
	BASE_bin = args.get("BASE_bin")

	base_par = []
	UGA = [i for i in range(DIM_H)]

	j=-1
	p=-1

	for i in range(DIM_H):

		if UGA[i] == "hola" :
			base_par.append((-1,-1,-1,-1))
			continue

		p_state,index_p_state = parity(BASE_bin[i],**args)

		UGA[index_p_state] 	  = "hola"
		
		j += 1

		if i == index_p_state:
			base_par.append((j,i,index_p_state,-1))
		else:
			p+= 1
			base_par.append((j,i,index_p_state,p))

	j += 1
	print(j)

	return base_par, j


def bose_Hamiltonian_parity(H_tmp,**args):

	mat_type    = args.get("mat_type")
	b_p_inp	    = args.get("parity_index")
	DIM_H 	    = np.int(args.get("DIM_H"))
	sim_sec_len = args.get("sim_sec_len")

	if mat_type == 'Sparse':
		H_tmp = H_tmp.todense()

	b_p   = np.asarray(b_p_inp)
	DX    = sim_sec_len

	H_par = np.matlib.zeros((DIM_H,DIM_H), dtype=np.float)

	for i in range(len(b_p)):

		if b_p[i,0] < 0:
			continue

		ind = b_p[i,0]

		if b_p[i,1] == b_p[i,2]:
			H_par[ind,:] 	+= H_tmp[int(b_p[i,1]),:]
		else:
			H_par[ind,:]    += (1/np.sqrt(2))*H_tmp[int(b_p[i,1]),:] 
			H_par[ind,:]    += (1/np.sqrt(2))*H_tmp[int(b_p[i,2]),:]
			H_par[DX,:]     += (1/np.sqrt(2))*H_tmp[int(b_p[i,1]),:]
			H_par[DX,:]     -= (1/np.sqrt(2))*H_tmp[int(b_p[i,2]),:]
			DX += 1

	H_dense = H_par
	H_par = np.matlib.zeros((DIM_H,DIM_H), dtype=np.float)

	DX = sim_sec_len

	for i in range(len(b_p)):
		
		if b_p[i,0] < 0:
			continue

		ind = b_p[i,0]

		if b_p[i,1] == b_p[i,2]:		
			H_par[:,ind]  += H_dense[:,int(b_p[i,1])]
		else:			
			H_par[:,ind]  += (1/np.sqrt(2))*H_dense[:,int(b_p[i,1])] 
			H_par[:,ind]  += (1/np.sqrt(2))*H_dense[:,int(b_p[i,2])]
			H_par[:,DX]   += (1/np.sqrt(2))*H_dense[:,int(b_p[i,1])] 
			H_par[:,DX]   -= (1/np.sqrt(2))*H_dense[:,int(b_p[i,2])]
			DX += 1

	if mat_type == 'Sparse':
		H_par = csc_matrix(H_par, shape=(DIM_H,DIM_H), dtype=np.double)

	return H_par


def vectors_parity_symmetrize(V1,**args):

	b_p_inp	 = args.get("parity_index")
	DIM_H 	 = np.int(args.get("DIM_H"))
	
	V0=np.transpose(V1)

	V   = np.matlib.zeros(np.shape(V0), dtype=np.float)
	
	b_p = np.asarray(b_p_inp)	
	DX  = args.get("sim_sec_len")


	for i in range(len(b_p)):

		if b_p[i,0] < 0:
			continue

		if b_p[i,1] == b_p[i,2]:
			V[:,int(b_p[i,1])] += V0[:,i]

		else:
			
			V[:,int(b_p[i,1])]  += np.sqrt(2)/2*V0[:,i]
			V[:,int(b_p[i,2])]  += np.sqrt(2)/2*V0[:,i]
			
			V[:,int(b_p[i,1])]  += np.sqrt(2)/2*V0[:,DX]
			V[:,int(b_p[i,2])]  -= np.sqrt(2)/2*V0[:,DX]

			DX += 1

	Vf = np.asarray(np.transpose(V))

	return Vf


def bose_Hamiltonian_parity_fast(**args):

	DIM_H 	  = np.int(args.get("DIM_H"))
	BASE_bin  = args.get("BASE_bin")
	BASE_bose = args.get("BASE_bose")
	mat_type  = args.get("mat_type")
	b_p_inp	  = args.get("parity_index")
	b_p       = np.asarray(b_p_inp)

	len_sym   = args.get("sim_sec_len")
	len_b_p   = len(b_p)

	len_asym  = DIM_H - len_sym

	X0_s = []
	Y0_s = []
	A0_s = []

	X0_a = []
	Y0_a = []
	A0_a = []



	for i in range(len_b_p):

		if b_p[i,0] < 0:
			continue

		X_1,Y_1,A_1 = ham.evaluate_ham( b_p[i,1] ,**args)
		X_2,Y_2,A_2 = ham.evaluate_ham( b_p[i,2] ,**args)

		X = [item for sublist in [X_1,X_2] for item in sublist]
		Y = [item for sublist in [Y_1,Y_2] for item in sublist]
		A = [item for sublist in [A_1,A_2] for item in sublist]

		for j in range(len(A)):

			state_X_0    = BASE_bin[X[j]]
			state_Y_0    = BASE_bin[Y[j]]

			state_X_rev  = state_X_0[::-1]
			state_Y_rev  = state_Y_0[::-1]


			ind_X 		 = X[j]
			ind_X_rev    = ff.get_index(state_X_rev,**args)
			
			ind_Y 		 = Y[j]
			ind_Y_rev    = ff.get_index(state_Y_rev,**args)	


			ind_col_X	 = b_p[min(ind_X,ind_X_rev),0]
			ind_col_Y	 = b_p[min(ind_Y,ind_Y_rev),0]
		
				
			coef_s = 2

		##.... SYM SEC

			if ind_X == ind_X_rev:

				ind_col_X = b_p[ind_X,0]
				coef_s 	  = 2*np.sqrt(2)

			X0_s.append(ind_col_X)	

			if ind_Y == ind_Y_rev:

				ind_col_Y = b_p[ind_Y,0]
				coef_s 	  = 2/np.sqrt(2)

				if ind_X == ind_X_rev:		

					coef_s = 2

			Y0_s.append(ind_col_Y) 
			A0_s.append(A[j]/coef_s)

		##.... A_SYM SEC
		
		for j in range(len(A_1)):

			state_X_0    = BASE_bin[X[j]]
			state_Y_0    = BASE_bin[Y[j]]

			state_X_rev  = state_X_0[::-1]
			state_Y_rev  = state_Y_0[::-1]


			ind_X 		 = X[j]
			ind_X_rev    = ff.get_index(state_X_rev,**args)
			
			ind_Y 		 = Y[j]
			ind_Y_rev    = ff.get_index(state_Y_rev,**args)	


			ind_col_X	 = b_p[min(ind_X,ind_X_rev),3]
			ind_col_Y	 = b_p[min(ind_Y,ind_Y_rev),3]

			if Y[j] > ind_Y_rev:	
				coef_a = -1

			elif Y[j] < ind_Y_rev:
				coef_a = 1
			
			else:
				continue

			X0_a.append(ind_col_X+len_sym)					
			Y0_a.append(ind_col_Y+len_sym) 
			A0_a.append(A[j]*coef_a)


	#X = [item for sublist in [X0_a,X0_s] for item in sublist]
	#Y = [item for sublist in [Y0_a,Y0_s] for item in sublist]
	#A = [item for sublist in [A0_a,A0_s] for item in sublist]

	#Hamiltonian = csc_matrix((A, (X,Y)), shape=(DIM_H,DIM_H), dtype=np.double)

	Hamiltonian_sym  = csc_matrix((A0_s, (X0_s,Y0_s)), shape=(len_sym,len_sym), dtype=np.double)
	Hamiltonian_asym = csc_matrix((A0_a, (X0_a,Y0_a)), shape=(len_asym,len_asym), dtype=np.double)


	if mat_type == 'Dense':

		#Hamiltonian  = csc_matrix.todense(Hamiltonian)

		Hamiltonian_sym  = csc_matrix.todense(Hamiltonian_sym)
		Hamiltonian_asym = csc_matrix.todense(Hamiltonian_asym)

	#ff.print_matrix(Hamiltonian)

	#ff.print_matrix(Hamiltonian_sym)
	#ff.print_matrix(Hamiltonian_asym)

	return Hamiltonian_sym, Hamiltonian_asym
