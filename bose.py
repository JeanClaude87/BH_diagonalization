import os
import profile
import numpy as np
import scipy as sp
from scipy.sparse import csc_matrix
from scipy.sparse import linalg
from numpy import linalg as LA
import time


import hamiltonian        as ham
import hamiltonian_parity as ham_par
import function           as ff
import observables        as ob
import time_evolution	  as t_ev


np.set_printoptions(precision=3)

ll_inp 			 = 9
nn_inp 			 = 9
BC_inp 			 = 0			# 0 is periodic
t_inp  			 = -1
U_inp  			 = -5
mat_type_inp     = 'Sparse' 	#'Sparse' #.... deafault Dense
parity_inp       = 'True'		#.... deafault False
n_diag_state_inp = 1
cores_num_inp    = 1


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
	"LOCAL" : os.path.abspath('.'),
	"parity"   : parity_inp,
	}

n_diag_state 	= Constants_dictionary.get("n_diag_state")

Constants_dictionary["tab_fact"]     = ff.fact_creation(**Constants_dictionary)

DIM_H 			= ff.hilb_dim(nn_inp, ll_inp, Constants_dictionary.get("tab_fact"))
Constants_dictionary["DIM_H"]        = DIM_H 
Constants_dictionary["hilb_dim_tab"] = ff.hilb_dim_tab(**Constants_dictionary)

print('Hilbert space Dimension:', Constants_dictionary.get("DIM_H"))
print('ll', ll_inp, 'nn', nn_inp)

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

quit()

E,V   = ham.diagonalization(Hamiltonian, **Global_dictionary)

















dt       = 0.1
step_num = 2500
t_i 	 = 0
t_f 	 = dt*step_num

#'''
part_ind = [4,4] 		# say in which site you want a particle
psi_4 = t_ev.inital_state(part_ind, **Global_dictionary)

part_ind = [9,9]
psi_9 = t_ev.inital_state(part_ind, **Global_dictionary)

psi_0 = 1/np.sqrt(2)*psi_4+1/np.sqrt(2)*psi_9
#'''

#part_ind = [4,4,9,9] 		# say in which site you want a particle
#psi_0 = t_ev.inital_state(part_ind, **Global_dictionary)

A        = -1.0J*Hamiltonian

psit     = linalg.expm_multiply(A, psi_0, start=t_i, stop=t_f, num=step_num+1, endpoint=True)
psit_par = ham_par.vectors_parity_symmetrize( psit.T, **Global_dictionary)

#ob.Export_Observable_time(psit_par,dt,'2+2_dens_t.dat',**Global_dictionary)

ob.Export_Fidelity(psit_par,dt,'fidelity.dat',**Global_dictionary)




quit()


part_ind = [4,4,4,4,9,9,9,9] 		# say in which site you want a particle
psi_0 = t_ev.inital_state(part_ind, **Global_dictionary)

psit     = linalg.expm_multiply(A, psi_0, start=t_i, stop=t_f, num=step_num+1, endpoint=True)
psit_par = ham_par.vectors_parity_symmetrize( psit.T, **Global_dictionary)

ob.Export_Observable_time(psit_par,dt,'22_dens_t.dat',**Global_dictionary)




