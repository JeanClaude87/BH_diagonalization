import os
import profile
import numpy as np
import time

import hamiltonian as ham
import function as ff
import observables as ob

np.set_printoptions(precision=5)

t1 = time.time()

ll_inp = 30
nn_inp = 3
BC_inp = 0
t_inp  = -1
U_inp  = -1
mat_type_inp = 'Sparse'

cores_num = 1


######............PREPARATION OF DICTIONARSS

Constants_dictionary = {}
Global_dictionary    = {}

Constants_dictionary = {
	"ll" : ll_inp, 
	"nn" : nn_inp,
	"BC" : BC_inp, 
	"t"  : t_inp ,
	"U"  : U_inp ,
	"cores_num" : cores_num,
	"mat_type" : mat_type_inp,
	"PATH_now" : os.path.abspath('.'),
	}

Constants_dictionary["tab_fact"] = ff.fact_creation(**Constants_dictionary)
Constants_dictionary["DIM_H"]    = ff.hilb_dim(nn_inp, ll_inp, Constants_dictionary.get("tab_fact"))

print('dimH', Constants_dictionary.get("DIM_H"))

Global_dictionary.update(Constants_dictionary)


BASE_bin, BASE_bose, CONF_tab = ff.Base_prep(**Constants_dictionary)

Global_dictionary["BASE_bin"]    = BASE_bin		#.......11100000, str
Global_dictionary["BASE_bose"]   = BASE_bose	#.......[3 0 0 0 0 0], numpy.ndarray
Global_dictionary["CONF_tab"]    = CONF_tab		#.......224, int

HOP_list     = ff.Hop_prep(**Constants_dictionary)

Global_dictionary["HOP_list"]  = HOP_list


Hamiltonian  = ham.bose_Hamiltonian(**Global_dictionary)


t2 = time.time()
print('Dt 1', t2-t1)

n_diag_state = 1

E,V = ham.diagonalization(Hamiltonian,n_diag_state,**Constants_dictionary)


t3 = time.time()
print('Dt 2', t3-t2)

for i in range(n_diag_state):

	dens   = ob.density( V[:,i],       **Global_dictionary)
	CdiCj  = ob.CdiCj(   V[:,i], dens, **Global_dictionary)

t4 = time.time()
print('Dt 3', t4-t3)

print('tot time', t4-t1)





