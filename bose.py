import sys
import os
import profile
import numpy as np
import scipy as sp
from scipy.optimize import curve_fit
from scipy.sparse import csc_matrix
from scipy.sparse import linalg as linalgS
from numpy import linalg as lin
import time
import matplotlib.pyplot as plt
from mpi4py import MPI

import hamiltonian        as ham
import hamiltonian_parity as ham_par
import function           as ff
import observables        as ob
import time_evolution	  as t_ev
import Hamiltonian_MPI	  as ham_MPI

np.set_printoptions(precision=6,suppress=True)

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

for nn_inp in [2]:#,3,4,5,6]:

	for ll_inp in [5]:

		if nn_inp == 2:	U_in = 5.0
		if nn_inp == 3:	U_in = 3.0	
		if nn_inp == 4:	U_in = 2.0	
		if nn_inp == 5:	U_in = 1.5	
		if nn_inp == 6:	U_in = 1.0						

#			for U_in in [0.5,1.0,3.0]:

		U_inp = -1.0*U_in
		
	for bar_inp in [0.05, 0.03, 0.01, 0.007, 0.005, 0.003, 0.001, 0.0007, 0.0005, 0.0003, 0.0001]:
	#for bar_inp in np.arange(0.,0.5,0.01)::
	#for bar_inp in [0.003]:	

		if nn_inp == 2:	bar_inp = 0.007
		if nn_inp == 3:	bar_inp = 0.003
		if nn_inp == 4:	bar_inp = 0.001	
		if nn_inp == 5:	bar_inp = 0.001	
		if nn_inp == 6:	bar_inp = 0.0007	

		#for flux_inp in np.arange(0.,0.5,0.01):

		flux_inp 		= 0.0
		flux_inp_1 		= 1.0
		flux_inp_psi0 	= 0.0
		flux_inp_t 		= 0.5		
		BC_inp 			= 0			# 0 is periodic

		#mat_type_inp     = 'Dense' 	#'Sparse' #.... default Dense
		mat_type_inp     = 'Sparse' 	#.... default Dense
		parity_inp       = 'False'		#.... default False
		n_diag_state_inp = 1 #ll_inp+4
		cores_num_inp    = 2
		t_inp  			 = -1*np.exp(-2*np.pi*1j*flux_inp/ll_inp)

		dt 		 = 10
		step_num = 1000

		#t max 4000

		if mat_type_inp == None:
			mat_type_inp = 'Sparse'

		######............PREPARATION OF DICTIONARSS

		Constants_dictionary = {}
		Global_dictionary    = {}

		zero = 0.0

		Constants_dictionary = {
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
			"step_num"	: step_num
			}

		n_diag_state 	= Constants_dictionary.get("n_diag_state")

		Constants_dictionary["tab_fact"]     = ff.fact_creation(**Constants_dictionary)

		DIM_H 			= ff.hilb_dim(nn_inp, ll_inp, Constants_dictionary.get("tab_fact"))

		Constants_dictionary["DIM_H"]        = DIM_H 
		Constants_dictionary["hilb_dim_tab"] = ff.hilb_dim_tab(**Constants_dictionary)

		if COMM.rank == 0:
			print('Hilbert space Dimension:', Constants_dictionary.get("DIM_H"))
			print('ll', ll_inp, 'nn', nn_inp)

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

		CDC 		 = ob.CdiCj_creation(**Global_dictionary)		

#		print(ob.C_d(1,**Global_dictionary))

#		profile.run('CDCDCC 		 = ob.CdCdCC_creation(**Global_dictionary)')			
#		print("mat", 'll', ll_inp, 'nn', nn_inp)

		Global_dictionary["CDC_matrix"]   = CDC
#		Global_dictionary["CDCDCC_matrix"]   = CDCDCC





#################### HAMILTONIAN 1  CREATION OMEGA = 0

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


				Hamiltonian_1 = csc_matrix((A1, (X1,Y1)), shape=(DIM_H,DIM_H), dtype=np.complex)

				if mat_type == 'Dense':

					Hamiltonian_1 = csc_matrix.todense(Hamiltonian_1)

			print("B0", 'fl', flux_inp, 'ba', Global_dictionary.get("bar"), 'U', Global_dictionary.get("U"))

			E0,V_cat_0  = ham.diagonalization(Hamiltonian_1, **Global_dictionary)



#################### HAMILTONIAN 2 CREATION OMEGA = 1.0


		COMM.Barrier()	

		t_inp_1  	= -1.0*np.exp(-2*np.pi*1j*flux_inp_1/ll_inp)
		
		Constants_dictionary = { 
		"t"  : t_inp_1
		}
		
		Global_dictionary.update(Constants_dictionary)
		
		COMM.Barrier()

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


			Hamiltonian_2 = csc_matrix((A1, (X1,Y1)), shape=(DIM_H,DIM_H), dtype=np.complex)
			
			if mat_type == 'Dense':

				Hamiltonian_2 = csc_matrix.todense(Hamiltonian_2)

			print("B1", 'fl', flux_inp_1, 'ba', Global_dictionary.get("bar"), 'U', Global_dictionary.get("U"))

			E1,V_cat_1  = ham.diagonalization(Hamiltonian_2, **Global_dictionary)

	
#################### HAMILTONIAN 3 CREATION OMEGA = 0.0 con barriera per psi0


		COMM.Barrier()	

		t_inp_t_psi0  	= -1.0*np.exp(-2*np.pi*1j*flux_inp_psi0/ll_inp)

		
		Constants_dictionary = { 
		"t"  : t_inp_t_psi0,
		"bar": bar_inp
		}
		
		Global_dictionary.update(Constants_dictionary)
		COMM.Barrier()

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

			Hamiltonian_3 = csc_matrix((A1, (X1,Y1)), shape=(DIM_H,DIM_H), dtype=np.complex)
			
			if mat_type == 'Dense':

				Hamiltonian_3 = csc_matrix.todense(Hamiltonian_3)
			
			print("H0", 'fl', flux_inp_psi0, 'ba', Global_dictionary.get("bar"), 'U', Global_dictionary.get("U"))
		
			Egs,V0  = ham.diagonalization(Hamiltonian_3, **Global_dictionary)
			

#################### HAMILTONIAN EVOLUTION CREATION OMEGA = 0.5

		COMM.Barrier()

		t_inp_t  	= -1.0*np.exp(-2*np.pi*1j*flux_inp_t/ll_inp)

		if COMM.rank == 0:
			LOCAL = os.path.abspath('.')
			directory = LOCAL+os.sep #+str('DATA')+os.sep+str('N_')+str(nn_inp)+os.sep+str('L_')+str(ll_inp)+os.sep+str('U_')+str(U_inp)+os.sep+str('bar_')+"%.9f"%bar_inp+os.sep+str('fl_')+str(flux_inp)+os.sep+str('fl-tt_')+str(flux_inp_t)+os.sep

			if not os.path.exists(directory):
				os.makedirs(directory)


		mat_type_inp     = 'Dense' 	#.... default Dense

		#if mat_type_inp == 'Dense':
		#	n_diag_state_inp = DIM_H
		
		Constants_dictionary = { 
		"n_diag_state" : n_diag_state_inp,
		"t"  : t_inp_t,
		"bar": bar_inp,
		"mat_type" 	: mat_type_inp
		}

		Global_dictionary.update(Constants_dictionary)
		
		COMM.Barrier()

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


			Hamiltonian_ev = csc_matrix((A1, (X1,Y1)), shape=(DIM_H,DIM_H), dtype=np.complex)

			if mat_type == 'Dense':

				Hamiltonian_ev = csc_matrix.todense(Hamiltonian_ev)

			print("ev", 'fl', flux_inp_t, 'ba', Global_dictionary.get("bar"), 'U', Global_dictionary.get("U"))


####################	TIME EVOLUTION -->>

		if COMM.rank == 0:
		
			psi_0 = V0
			psit = t_ev.time_evolution(psi_0, Hamiltonian_ev, **Global_dictionary)


####################	OBSERVABLES -->> 
			

			directory = os.sep+'dati'+os.sep+'L_'+str(ll_inp)+os.sep+'N_'+str(nn_inp)+os.sep+'U_'+str(U_inp)+os.sep+'bb_'+str(bar_inp)

			Dstep = 1

	#		CCDD    = ob.CdCdCC_t  (psit, Dstep, **Global_dictionary)
	#		CD      = ob.CdiCj_t   (psit, Dstep, **Global_dictionary)
			current = ob.corrente_t(psit, Dstep, **Global_dictionary)

	#		ob.Export_Observable(CCDD, directory,    'densdens.dat', **Global_dictionary)
	#		ob.Export_Observable(CD,   directory,    'dens.dat',     **Global_dictionary)
			ob.Export_Observable(current, directory, 'corrente.dat', **Global_dictionary)								
			
			ob.Export_Fidelity_CAT_s(psit, V_cat_0, V_cat_1, directory, 'fidelity_cat_s.dat',**Global_dictionary)
			ob.Export_Fidelity_CAT_a(psit, V_cat_0, V_cat_1, directory, 'fidelity_cat_a.dat',**Global_dictionary)			
			ob.Export_Fidelity(psit, V_cat_0,   directory, 'fidelity_0.dat',**Global_dictionary)
			ob.Export_Fidelity(psit, V_cat_1,   directory, 'fidelity_1.dat',**Global_dictionary)


quit()


