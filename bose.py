import os

import numpy as np
from math import factorial
import math
import itertools
from scipy.sparse import csc_matrix
from scipy.sparse import linalg
from numpy import linalg as LA
import time

import hamiltonian as ham
import function as ff
import observables as ob

np.set_printoptions(precision=2)

t1 = time.clock()

ll_inp = 14
nn_inp = 5
BC_inp = 0
t_inp  = -1
U_inp  = -1
mat_type_inp = 'Sparse'


######............PREPARATION OF DICTIONARSS

Constants_dictionary = {}
Global_dictionary    = {}

Constants_dictionary = {
	"ll" : ll_inp, 
	"nn" : nn_inp,
	"BC" : BC_inp, 
	"t"  : t_inp ,
	"U"  : U_inp ,
	"mat_type" : mat_type_inp,
	"PATH_now" : os.path.abspath('.'),
	}

Constants_dictionary["tab_fact"] = ff.fact_creation(**Constants_dictionary)
Constants_dictionary["DIM_H"]    = ff.hilb_dim(nn_inp, ll_inp, Constants_dictionary.get("tab_fact"))

print(Constants_dictionary.get("DIM_H"))

Global_dictionary.update(Constants_dictionary)


BASE_bin, BASE_bose, CONF_tab = ff.Base_prep(**Constants_dictionary)

Global_dictionary["BASE_bin"]  = BASE_bin
Global_dictionary["BASE_bose"] = BASE_bose
Global_dictionary["CONF_tab"]  = CONF_tab

HOP_list     = ff.Hop_prep(**Constants_dictionary)

Global_dictionary["HOP_list"]  = HOP_list


Hamiltonian  = ham.bose_Hamiltonian(**Global_dictionary)

t2 = time.clock()
print('Dt 1', t2-t1)

n_diag_state = 3

E,V = ham.diagonalization(Hamiltonian,n_diag_state,**Constants_dictionary)

for i in range(n_diag_state):

	dens0  = ob.density(   V[:,i],BASE_bose)
	#print(dens0)

t3 = time.clock()

print('Dt 2', t3-t2)






