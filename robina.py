import sys
import os
import profile
import numpy as np
import scipy as sp
from scipy.optimize import curve_fit
from scipy.sparse import csc_matrix
from scipy.sparse import linalg
from numpy import linalg as LA
import time
import matplotlib.pyplot as plt
from mpi4py import MPI

import hamiltonian        as ham
import hamiltonian_parity as ham_par
import function           as ff
import observables        as ob
import time_evolution	  as t_ev
import Hamiltonian_MPI	  as ham_MPI

np.set_printoptions(precision=8,suppress=True)

COMM = MPI.COMM_WORLD

if COMM.rank == 0:
	t1 = time.time()


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

for ll_inp in [30]:

	nn_inp 	=	4
	#ll_inp 	=	12
	bar_inp =	0.00


	BC_inp 			 = 0			# 0 is periodic

	mat_type_inp     = 'Sparse'#Dense' 	#.... default Dense
	parity_inp       = 'False'		#.... default False
	n_diag_state_inp = 1
	cores_num_inp    = 2

		# term    coeff	t*exp(+2*pi*iu*omegai/Sites)

	if mat_type_inp == None:
		mat_type_inp = 'Sparse'

	######............PREPARATION OF DICTIONARSS

	Constants_dictionary = {}
	Global_dictionary    = {}

	Constants_dictionary = {
		"ll" : ll_inp, 
		"nn" : nn_inp,
		"BC" : BC_inp, 
		"bar": bar_inp ,
		"n_diag_state"  : n_diag_state_inp,
		"cores_num" : cores_num_inp,
		"mat_type" 	: mat_type_inp,
		"LOCAL" 	: os.path.abspath('.'),
		"parity"   	: parity_inp,
		}

	n_diag_state 	= Constants_dictionary.get("n_diag_state")

	Constants_dictionary["tab_fact"]     = ff.fact_creation(**Constants_dictionary)

	DIM_H 			= ff.hilb_dim(nn_inp, ll_inp, Constants_dictionary.get("tab_fact"))

	Constants_dictionary["DIM_H"]        = DIM_H 
	Constants_dictionary["hilb_dim_tab"] = ff.hilb_dim_tab(**Constants_dictionary)

	#if COMM.rank == 0:
	#	print('Hilbert space Dimension:', Constants_dictionary.get("DIM_H"))
	#	print('ll', ll_inp, 'nn', nn_inp)

	Global_dictionary.update(Constants_dictionary)

	COMM.Barrier()
	#print(Global_dictionary.get("t"))

	if COMM.rank == 0:

		BASE_bin, BASE_bose, CONF_tab = ff.Base_prep(**Constants_dictionary)

		Global_dictionary["BASE_bin"]    = BASE_bin		#.......11100000, str
		Global_dictionary["BASE_bose"]   = BASE_bose	#.......[3 0 0 0 0 0], numpy.ndarray
		Global_dictionary["CONF_tab"]    = CONF_tab		#.......224, int
		#print(BASE_bose)

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


	HOP_list     = ff.Hop_prep(**Constants_dictionary)

	Global_dictionary["HOP_list"]  = HOP_list

#	CDC 			= ob.CdiCj_creation(**Global_dictionary)
#	Global_dictionary["CDC_matrix"]   = CDC


	for U_0 in [0.5]: #np.arange(0.,0.5,0.01):

		U_inp = -1.0*U_0

		for flux_inp in np.arange(0.,0.52,0.02):
			
			t_inp = -1*np.exp(-2*np.pi*1j*flux_inp/ll_inp)
			
			Constants_dictionary = {
				"U"  : U_inp ,
				"flux": flux_inp,
				"t"  : t_inp
				}

			Global_dictionary.update(Constants_dictionary)

			#################### HAMILTONIAN 1  CREATION OMEGA = 0

			#print("V_cat_0", flux_inp , Global_dictionary.get("U"), Global_dictionary.get("bar"))

			if Constants_dictionary.get("parity") == 'True':

				Global_dictionary["parity_index"], Constants_dictionary["sim_sec_len"] = ham_par.base_parity(**Global_dictionary)
				Global_dictionary.update(Constants_dictionary)

				#print('I do parity!! ')
				if COMM.rank == 0:
				
					Hamiltonian = ham_par.bose_Hamiltonian_parity_fast(**Global_dictionary)

			else:

				if COMM.rank == 0:

					mat_type = Global_dictionary.get("mat_type")
					jobs = list(range(DIM_H))
					jobs = ham_MPI.split(jobs, COMM.size)

				else:
					jobs = None

				COMM.Barrier()

				jobs = COMM.scatter(jobs, root=0)

				XX = []
				YY = []
				AA = []

				for i in jobs:
					res = ham.evaluate_ham(i, **Global_dictionary)

					XX.append(res[0])
					YY.append(res[1])
					AA.append(res[2])

				COMM.Barrier()

				XX0 = MPI.COMM_WORLD.gather( XX, root=0)
				YY0 = MPI.COMM_WORLD.gather( YY, root=0)
				AA0 = MPI.COMM_WORLD.gather( AA, root=0)

				COMM.Barrier()

				if COMM.rank == 0:

					X0 = [item for sublist in XX0 for item in sublist]
					Y0 = [item for sublist in YY0 for item in sublist]
					A0 = [item for sublist in AA0 for item in sublist]

					X1 = [item for sublist in X0 for item in sublist]
					Y1 = [item for sublist in Y0 for item in sublist]
					A1 = [item for sublist in A0 for item in sublist]


					Hamiltonian_0 = csc_matrix((A1, (X1,Y1)), shape=(DIM_H,DIM_H), dtype=np.complex)
					
					if mat_type == 'Dense':

						Hamiltonian = csc_matrix.todense(Hamiltonian_0)

			if COMM.rank == 0:

				E0,V_cat_0  = ham.diagonalization(Hamiltonian_0, **Global_dictionary)

				#3cc = ob.CdiCj(V_cat_0[:,0], **Global_dictionary)
				
				print(ll_inp, U_inp, flux_inp, E0[0])

quit()


