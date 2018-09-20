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

	for i in range(DIM_H):

		if UGA[i] == "hola" :
			continue

		p_state,index_p_state = parity(BASE_bin[i],**args)

		UGA[index_p_state] 	  = "hola"
		
		base_par.append((i,index_p_state))

	return base_par


def bose_Hamiltonian_parity(H_tmp,**args):

	mat_type = args.get("mat_type")
	b_p_inp	 = args.get("parity_index")
	DIM_H 	 = np.int(args.get("DIM_H"))

	if mat_type == 'Sparse':
		H_tmp = H_tmp.todense()


	b_p   = np.asarray(b_p_inp)
	DX    = len(b_p)

	H_par = np.matlib.zeros((DIM_H,DIM_H), dtype=np.float)

	for i in range(len(b_p)):

		if b_p[i,0] == b_p[i,1]:
			H_par[i,:] 	+= H_tmp[int(b_p[i,0]),:]
		else:
			H_par[i,:]  += (1/np.sqrt(2))*H_tmp[int(b_p[i,0]),:] 
			H_par[i,:]  += (1/np.sqrt(2))*H_tmp[int(b_p[i,1]),:]
			H_par[DX,:] += (1/np.sqrt(2))*H_tmp[int(b_p[i,0]),:]
			H_par[DX,:] -= (1/np.sqrt(2))*H_tmp[int(b_p[i,1]),:]
			DX += 1

	H_dense = H_par
	H_par = np.matlib.zeros((DIM_H,DIM_H), dtype=np.float)

	DX      = len(b_p)

	for i in range(len(b_p)):

		if b_p[i,0] == b_p[i,1]:
			H_par[:,i]  += H_dense[:,int(b_p[i,0])]
		else:
			H_par[:,i]  += (1/np.sqrt(2))*H_dense[:,int(b_p[i,0])] 
			H_par[:,i]  += (1/np.sqrt(2))*H_dense[:,int(b_p[i,1])]
			H_par[:,DX] += (1/np.sqrt(2))*H_dense[:,int(b_p[i,0])] 
			H_par[:,DX] -= (1/np.sqrt(2))*H_dense[:,int(b_p[i,1])]
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
	DX  = len(b_p)


	for i in range(len(b_p)):

		if b_p[i,0] == b_p[i,1]:
			V[:,int(b_p[i,0])] += V0[:,i]

		else:
			
			V[:,int(b_p[i,0])]  += np.sqrt(2)/2*V0[:,i]
			V[:,int(b_p[i,1])]  += np.sqrt(2)/2*V0[:,i]
			
			V[:,int(b_p[i,0])]  += np.sqrt(2)/2*V0[:,DX]
			V[:,int(b_p[i,1])]  -= np.sqrt(2)/2*V0[:,DX]

			DX += 1

	Vf = np.asarray(np.transpose(V))

	return Vf


def bose_Hamiltonian_parity_fast(**args):

	DIM_H 	 = np.int(args.get("DIM_H"))
	BASE_bin = args.get("BASE_bin")
	mat_type = args.get("mat_type")
	b_p_inp	 = args.get("parity_index")
	b_p = np.asarray(b_p_inp)

	len_sym  = len(b_p)
	len_asym = DIM_H - len_sym

#	H_par_sym  = np.matlib.zeros((len_sym,len_sym), dtype=np.float)
#	H_par_asym = np.matlib.zeros((len_asym,len_asym), dtype=np.float)


	X0_s = []
	Y0_s = []
	A0_s = []


	X0_a = []
	Y0_a = []
	A0_a = []

	for i in range(len_sym):

		if b_p[i,0] == b_p[i,1]:

			state     = BASE_bin[b_p[i,0]]

			X,Y,A = ham.evaluate_ham( b_p[i,0] ,**args)

			for j in range(len(A)):
				
				state_0    = BASE_bin[Y[j]]
				state_rev  = state_0[::-1]

				ind_0   = ff.get_index(state_0,**args)
				ind_rev = ff.get_index(state_rev,**args)
	

				print(b_p[i,0],ind_0,ind_rev,min(ind_0,ind_rev))

#				print(X[i],Y[i],A[i])
	
#			X0_s.append(A)
#			Y0_s.append(B)
#			A0_s.append(C)

		else:
			state0 = BASE_bin[b_p[i,0]]
			state1 = BASE_bin[b_p[i,0]]


#	X1_s = [item for sublist in X0_s for item in sublist]
#	Y1_s = [item for sublist in Y0_s for item in sublist]
#	A1_s = [item for sublist in A0_s for item in sublist]

#	Hamiltonian_sym = csc_matrix((A1_s, (X1_s,Y1_s)), shape=(len_sym,len_sym), dtype=np.double)

#	if mat_type == 'Dense':

#		Hamiltonian_sym = csc_matrix.todense(Hamiltonian)

#	ff.print_matrix(Hamiltonian)






	return 0

