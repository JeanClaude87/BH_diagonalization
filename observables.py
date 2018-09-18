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


## .................................................................
## ....................OBSERVABLES..................................
## .................................................................

#def density():

#..................................................dens
def density(V,BASE_bose):

	den   = np.dot(np.transpose(V**2),BASE_bose)

	return den

def OUTER_creation(A):
	
	Dim = len(A)
	L = len(A[0])
	B = np.zeros((Dim,L,L), dtype=np.float)

	for i in range(Dim):
		B[i] = np.outer(A[i],A[i])
	return B

def NiNj(V,CORR_BASE):

	NN = np.einsum('n,nij -> ij', V**2, CORR_BASE)

	return NN

def NfixNr(V,i,CORR_BASE):

	NiN = np.einsum('n,nj -> j', V**2, CORR_BASE[:,i])

	return NiN

def CdiCj(**args):

	state    = np.asarray(args.get("BASE_bose"))
	DIM_H 	 = args.get("DIM_H")
	ll  	 = args.get("ll")


	for i, j in itertools.product(range(ll), range(ll)):

		t1=time.time()

		BB = np.zeros((DIM_H,DIM_H))

		operator    = [0 for i in range(ll)]
		operator[i] = +1
		operator[j] = -1

		for psi0 in range(DIM_H):

			XX = []
			YY = []
			AA = []

			a = state[psi0] + operator

			for psi1 in range(DIM_H):

				if np.array_equal(state[psi1],a):
					
					XX.append(psi0)
					YY.append(psi1)
					AA.append(1)

					BB[psi0,psi1] = 1

					#print(i,'i', state[i],'.....',j,'j', state[j], '...', a)

		print(i,j)

	t2=time.time()
	
	return t2-t1









