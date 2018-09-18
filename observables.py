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

	state    = np.asarray(args.get("BASE_bin"))
	DIM_H 	 = args.get("DIM_H")
	ll  	 = args.get("ll")
	nn  	 = args.get("nn")	

	BB = np.zeros((ll,ll,DIM_H,DIM_H))

	t1=time.time()

	for i, j in itertools.product(range(ll), range(ll)):

		for psi0 in range(DIM_H):

			a = state[psi0]

	t2=time.time()
	
	return t2-t1









