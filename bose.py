import os
import profile
import numpy as np
from scipy.sparse import csc_matrix
import time
from mpi4py import MPI

import cProfile

import hamiltonian        as ham
import hamiltonian_parity as ham_par
import function           as ff
import observables        as ob
import time_evolution	  as t_ev
import Hamiltonian_MPI	  as ham_MPI

np.set_printoptions(precision=5,suppress=True)

COMM = MPI.COMM_WORLD


#  n 	 Uc
# 2.0	4.0	 0.0
# 3.0	1.83 0.05
# 4.0 	1.15 0.05
# 5.0 	0.82 0.05
# 6.0 	0.63 0.05
# 7.0 	0.48 0.08
# 8.0 	0.45 0.08
# 9.0 	0.40 0.08
# 10.0 	0.32 0.08


for nn_inp in [4]:
		
	bar_inp = 0.0
	ll_inp  = 16
	U_in 	= 1.0	
	
	U_inp   = -1.0*U_in

	#for ll_inp in np.arange(8,38,4):

	BC_inp 			= 0			# 0 is periodic

	#mat_type_inp     = 'Dense' 	#'Sparse' #.... default Dense
	mat_type_inp     = 'Sparse' 	#.... default Dense
	parity_inp       = 'False'		#.... default False
	n_diag_state_inp = 2
	cores_num_inp    = 2
	t_inp  			 = -1


	t_start  = 0
	dt 		 = 10
	step_num = 500	#100
	Dstep = 1

	#t max 4000

	if mat_type_inp == None:
		mat_type_inp = 'Sparse'

	######............PREPARATION OF DICTIONARSS

	zero = 0.0

	Global_dictionary = {
		"ll" 		: ll_inp, 
		"nn" 		: nn_inp,
		"BC" 		: BC_inp, 
		"t"  		: t_inp ,
		"U"  		: U_inp ,
		"bar"		: zero,
		"n_diag_state"  : n_diag_state_inp,
		"cores_num" : cores_num_inp,
		"mat_type" 	: mat_type_inp,
		"LOCAL" 	: os.path.abspath('.'),
		"parity"   	: parity_inp,
		"dt"      	: dt,
		"step_num"	: step_num,
		"t_start"	: t_start
		}

	n_diag_state 	= Global_dictionary.get("n_diag_state")

	Global_dictionary["tab_fact"]     = ff.fact_creation(ll_inp,nn_inp) #**Global_dictionary)

	DIM_H 			= ff.hilb_dim(nn_inp, ll_inp, Global_dictionary.get("tab_fact"))

	Global_dictionary["DIM_H"]        = DIM_H 
	Global_dictionary["hilb_dim_tab"] = ff.hilb_dim_tab(**Global_dictionary)

	if COMM.rank == 0:
		print('Hilbert space Dimension:', Global_dictionary.get("DIM_H"))
		print('ll', ll_inp, 'nn', nn_inp)

	COMM.Barrier()

	if COMM.rank == 0:

		BASE_bin, BASE_bose, CONF_tab = ff.Base_prep(**Global_dictionary)

		Global_dictionary["BASE_bin"]    = BASE_bin		#.......11100000, str
		Global_dictionary["BASE_bose"]   = BASE_bose	#.......[3 0 0 0 0 0], numpy.ndarray
		Global_dictionary["CONF_tab"]    = CONF_tab		#.......224, int

	else:

		BASE_bin 		= None
		BASE_bose 		= None
		CONF_tab		= None

	COMM.Barrier()
	BASE_bin 		= COMM.bcast(BASE_bin,	root=0)
	BASE_bose 		= COMM.bcast(BASE_bose,	root=0)
	CONF_tab		= COMM.bcast(CONF_tab,	root=0)

	Global_dictionary["BASE_bin"]    = BASE_bin		#.......11100000, str
	Global_dictionary["BASE_bose"]   = BASE_bose	#.......[3 0 0 0 0 0], numpy.ndarray
	Global_dictionary["CONF_tab"]    = CONF_tab		#.......224, int

	HOP_list     = ff.Hop_prep(**Global_dictionary)
	Global_dictionary["HOP_list"]  = HOP_list

	N 		 	= ob.N_creation(**Global_dictionary)	
	CDC 		= ob.CdiCj_creation(**Global_dictionary)	

	Global_dictionary["CDC_matrix"]   = CDC
	Global_dictionary["N_matrix"]     = N

	
	################################################################################################
	################################################################################################
	############ HAMILTONIAN CAT 0

	Hint 	= ob.int_op     (		**Global_dictionary)

	U_inp 	= -1.0

	for fi in np.arange(0,0.1,0.001):
		#for bar in ff.crange(0,0.5,0.01):

		HH   = U_inp/2*Hint + ob.kinetik_op (fi, **Global_dictionary)
		E,V  = ham.diagonalization(HH, **Global_dictionary)

		print(E)
		print(fi)

	quit()






















