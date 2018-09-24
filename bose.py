import os
import profile
import numpy as np
import scipy as sp
from scipy.sparse import csc_matrix
import time

import hamiltonian        as ham
import hamiltonian_parity as ham_par
import function           as ff
import observables        as ob


np.set_printoptions(precision=4)

ll_inp = 10
nn_inp = 6
BC_inp = 0
t_inp  = -1
U_inp  = -1
mat_type_inp = 'Sparse' #.... deafault Dense
parity_inp   = 'True'	#.... deafault False
n_diag_state_inp = 1

cores_num_inp = 1


######............PREPARATION OF DICTIONARSS

Constants_dictionary = {}
Global_dictionary    = {}

Constants_dictionary = {
	"ll" : ll_inp, 
	"nn" : nn_inp,
	"BC" : BC_inp, 
	"t"  : t_inp ,
	"U"  : U_inp ,
	"n_diag_state"  : n_diag_state_inp,
	"cores_num" : cores_num_inp,
	"mat_type" : mat_type_inp,
	"PATH_now" : os.path.abspath('.'),
	"parity"   : parity_inp,
	}

n_diag_state=Constants_dictionary.get("n_diag_state")


Constants_dictionary["tab_fact"]     = ff.fact_creation(**Constants_dictionary)
Constants_dictionary["DIM_H"]        = ff.hilb_dim(nn_inp, ll_inp, Constants_dictionary.get("tab_fact"))
Constants_dictionary["hilb_dim_tab"] = ff.hilb_dim_tab(**Constants_dictionary)

print('Hilbert space Dimension:', Constants_dictionary.get("DIM_H"))

Global_dictionary.update(Constants_dictionary)

BASE_bin, BASE_bose, CONF_tab = ff.Base_prep(**Constants_dictionary)

Global_dictionary["BASE_bin"]    = BASE_bin		#.......11100000, str
Global_dictionary["BASE_bose"]   = BASE_bose	#.......[3 0 0 0 0 0], numpy.ndarray
Global_dictionary["CONF_tab"]    = CONF_tab		#.......224, int

HOP_list     = ff.Hop_prep(**Constants_dictionary)

Global_dictionary["HOP_list"]  = HOP_list

if Constants_dictionary.get("parity") == 'True':

	Global_dictionary["parity_index"], Constants_dictionary["sim_sec_len"] = ham_par.base_parity(**Global_dictionary)
	Global_dictionary.update(Constants_dictionary)

	print('I do parity!! ')

	Hamiltonian = ham_par.bose_Hamiltonian_parity_fast(**Global_dictionary)

else:

	Hamiltonian = ham.bose_Hamiltonian(**Global_dictionary)


E,V = ham.diagonalization(Hamiltonian, **Global_dictionary)


for i in range(n_diag_state):

	dens   = ob.density( V[:,i],       **Global_dictionary)

	print(dens)





