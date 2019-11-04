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


for nn_inp in [2,3,4]:

		if nn_inp == 2: ll_inp = 40
		if nn_inp == 3: ll_inp = 30
		if nn_inp == 4: ll_inp = 20

		ciao=0

	#	'''

		if nn_inp == 1:	U_in = 5.0
		if nn_inp == 2:	U_in = 5.0
		if nn_inp == 3:	U_in = 3.0	
		if nn_inp == 4:	U_in = 2.0	
		if nn_inp == 5:	U_in = 1.5	
		if nn_inp == 6:	U_in = 1.0						

	#for U_in in [-6, -4, -2, 0, 3, 6]:

		U_inp = -1.0*U_in

		'''

		if nn_inp == 1:	bar_inp = 0.007
		if nn_inp == 2:	bar_inp = 0.0085
		if nn_inp == 3:	bar_inp = 0.003
		if nn_inp == 4:	bar_inp = 0.0025	
		if nn_inp == 5:	bar_inp = 0.001	
		if nn_inp == 6:	bar_inp = 0.0007			

		'''

		for bar_inp in np.arange(0.001,0.01,0.0005):#[0.0035, 0.0045]:#[0.007, 0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 0.00005, 0.0005, 0.005, 0.05, 0.5, 5, 50]:


			#print(bar_inp, ll_inp)

			flux_inp_0 		= 0.0
			flux_inp_1 		= 1.0
			flux_inp_psi0 	= 0.0
			flux_inp_t 		= 0.5		
			BC_inp 			= 0			# 0 is periodic

			#mat_type_inp     = 'Dense' 	#'Sparse' #.... default Dense
			mat_type_inp     = 'Sparse' 	#.... default Dense
			parity_inp       = 'False'		#.... default False
			n_diag_state_inp = 1
			cores_num_inp    = 2
			t_inp  			 = -1


			t_start  = 0
			dt 		 = 10
			step_num = 600	#100
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

			#if COMM.rank == 0:
				#print('Hilbert space Dimension:', Global_dictionary.get("DIM_H"))
				#print('ll', ll_inp, 'nn', nn_inp)

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
			NN 		    = ob.NNm1_creation(**Global_dictionary)	

			CDC 		= ob.CdiCj_creation(**Global_dictionary)	

			Global_dictionary["CDC_matrix"]   = CDC
			Global_dictionary["N_matrix"]     = N
			Global_dictionary["NN_matrix"]  = NN

			directory = os.sep+'dati'+os.sep+'L_'+str(ll_inp)+os.sep+'N_'+str(nn_inp)+os.sep+'U_'+str(U_inp)+os.sep+'bb_'+str(bar_inp)

			if COMM.rank == 0:
			
				Hint   = ob.int_op(**Global_dictionary)

				for om in [0.0]: #np.arange(0,1,0.02):

					Kin	   = ob.kinetik_op (om,  **Global_dictionary)
					cu_0   = ob.corrente_op(om,  **Global_dictionary)
					fl_0   = ob.fluct_op   (cu_0,**Global_dictionary)

					matrix_h = Kin + U_inp/2*Hint + bar_inp*ob.bar_0(0,**Global_dictionary)

					E,V0  = ham.diagonalization( matrix_h , **Global_dictionary)

					V   = V0.T[0]
					V_c = np.conjugate(V)					

					cu = np.real(V_c.dot(cu_0.dot(V)))
					fl = np.real(V_c.dot(fl_0.dot(V)))

					#print(nn_inp, ll_inp, U_inp, om, E[0], cu, fl)

					ob.Export_Observable([cu,fl], 	directory, 't=0.dat', **Global_dictionary)




		################################################################################################
		################################################################################################
		############ HAMILTONIAN CAT 0

				Hint   		= U_inp/2*ob.int_op(**Global_dictionary)
				
				Hkin_0 		= ob.kinetik_op (0,  **Global_dictionary)
				Hkin_1 		= ob.kinetik_op (1,  **Global_dictionary)
				Hkin_05 	= ob.kinetik_op (0.5,  **Global_dictionary)

				Hba_0   	= bar_inp*ob.bar_0		(0,**Global_dictionary)

		############ HAMILTONIAN CAT omega = 0 no barrier

				HH_1   = Hint + Hkin_0
				HH_1   = csc_matrix(HH_1, shape = (DIM_H,DIM_H))

				E1,V_cat_0  = ham.diagonalization(HH_1, **Global_dictionary)

		############ HAMILTONIAN CAT omega = 1 no barrier

				HH_2   = Hint + Hkin_1
				HH_2   = csc_matrix(HH_2, shape = (DIM_H,DIM_H))

				E2,V_cat_1  = ham.diagonalization(HH_2, **Global_dictionary)

		############ HAMILTONIAN t=0 omega = 0, YES barrier

				HH_3   = Hint + Hkin_0 + Hba_0
				HH_3   = csc_matrix(HH_3, shape = (DIM_H,DIM_H))

				E3,V3  = ham.diagonalization(HH_3, **Global_dictionary)

		############ HAMILTONIAN t evolution omega = 1/2, YES barrier

				HH_ev  = Hint + Hkin_05 + Hba_0
				#HH_ev  = csc_matrix(HH_ev, shape = (DIM_H,DIM_H))
                
		############ time_evolution

				psi_0 = V3
				psit = t_ev.time_evolution(psi_0, HH_ev, **Global_dictionary)

		####################	OBSERVABLES -->> 			
				ll = ll_inp

				cu_op 		= ob.corrente_op(0,    **Global_dictionary)
				fl_op   	= ob.fluct_op   (cu_op,**Global_dictionary)

				t_vec          = range(0,step_num,Dstep)
				corrente_array = np.real(np.array([[t*dt+t_start, np.conjugate(psit[t]).dot(cu_op.dot(psit[t]))] for t in t_vec]))
				fisherin_array = np.real(np.array([[t*dt+t_start, np.conjugate(psit[t]).dot(cu_op.dot(psit[t])), np.conjugate(psit[t]).dot(fl_op.dot(psit[t]))] for t in t_vec]))

				#print('fish', fisherin_array)				

				ob.Export_Observable(fisherin_array,   	directory, 'fish_t.dat',     **Global_dictionary)
				ob.Export_Observable(corrente_array, 	directory, 'corrente.dat', **Global_dictionary)										

				ob.Export_Fidelity_CAT_s(psit, V_cat_0, V_cat_1, directory, 'fidelity_cat_s.dat',**Global_dictionary)
				ob.Export_Fidelity_CAT_a(psit, V_cat_0, V_cat_1, directory, 'fidelity_cat_a.dat',**Global_dictionary)			
				ob.Export_Fidelity(psit, V_cat_0,   directory, 'fidelity_0.dat',**Global_dictionary)
				ob.Export_Fidelity(psit, V_cat_1,   directory, 'fidelity_1.dat',**Global_dictionary)
				ob.Export_Fidelity(psit, psi_0,   directory, 'fidelity_psi0.dat',**Global_dictionary)


#quit()

