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
def density(V,**args):
	
	BASE_bose = args.get("BASE_bose")

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


def CdiCj(V, dens, **args):

	states   = args.get("BASE_bose")
	DIM_H 	 = args.get("DIM_H")
	ll  	 = args.get("ll")
	nn  	 = args.get("nn")	

	Cd  = np.power(np.remainder(states+1,nn+1),1/4)
	C   = np.power(states,1/4)
	V_c = np.conj(V)

	CdiCj = np.zeros((ll,ll), dtype=np.float)

	for i in range(ll):
		for j in range(i+1,ll):

				Cd_C = Cd[:,i] * C[:,j]		# V
				C_Cd = C[:,i]  * Cd[:,j] 	# V*

				A = V   * Cd_C
				B = V_c * C_Cd
				
				CdiCj[i,j] = np.dot(A[A!=0],B[B!=0])

				

	CdiCj += CdiCj.T
	CdiCj += np.diag(dens)



	return CdiCj









