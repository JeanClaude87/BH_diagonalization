import sys
import os
import profile
import numpy as np
import scipy as sp
from scipy.sparse import csc_matrix
from scipy.sparse import linalg
from numpy import linalg as LA
import time

from mpi4py import MPI

import hamiltonian        as ham
import hamiltonian_parity as ham_par
import function           as ff
import observables        as ob
import time_evolution	  as t_ev
import Hamiltonian_MPI	  as ham_MPI

np.set_printoptions(precision=3)

COMM = MPI.COMM_WORLD

if COMM.rank == 0:
	t1 = time.time()


nn_inp   = int(sys.argv[1])
ll_inp   = int(sys.argv[2])
U_inp    = -1.0*float(sys.argv[3])
flux_inp = float(sys.argv[4])

#nn_inp 		= 2
#ll_inp 		= 10
#U_inp  		= -1
#flux_inp		= 0

BC_inp 			 = 0			# 0 is periodic

mat_type_inp     = 'Sparse' 	#.... deafault Dense
parity_inp       = 'False'		#.... deafault False
n_diag_state_inp = 1
cores_num_inp    = 2
t_inp  			 = -1*np.exp(-2*np.pi*1j*flux_inp/ll_inp)


#for nn_inp in [2, 3, 4, 5]:
#	for ll_inp in [5, 10, 15, 20, 25]:
#		for U_inp in np.arange(0.1,3,0.1):
#			for flux_inp in np.arange(0.00,0.5,0.02):

if mat_type_inp == None:
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

#if COMM.rank == 0:

#print('Hilbert space Dimension:', Constants_dictionary.get("DIM_H"))
	#print('ll', ll_inp, 'nn', nn_inp)

Global_dictionary.update(Constants_dictionary)

COMM.Barrier()


if COMM.rank == 0:

	BASE_bin, BASE_bose, CONF_tab = ff.Base_prep(**Constants_dictionary)

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


HOP_list     = ff.Hop_prep(**Constants_dictionary)

Global_dictionary["HOP_list"]  = HOP_list

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


		Hamiltonian = csc_matrix((A1, (X1,Y1)), shape=(DIM_H,DIM_H), dtype=np.complex)
		
		if mat_type == 'Dense':

			Hamiltonian = csc_matrix.todense(Hamiltonian)

if COMM.rank == 0:

	E,V   = ham.diagonalization(Hamiltonian, **Global_dictionary)

	#dens  = ob.density(V, **Global_dictionary)
	#CC    = ob.CdiCj  (V, dens, **Global_dictionary)
	
	nn_cor = ob.NiNj(V,**Global_dictionary)

	directory = 'DATA'+os.sep+'N_'+str(nn_inp)+os.sep+'L_'+str(ll_inp)+os.sep+'U_'+str(U_inp)+os.sep+'Om_'+str(flux_inp)
	LOCAL 	  = Constants_dictionary.get("LOCAL")
	PATH_now  = LOCAL+os.sep+directory+os.sep

	if not os.path.exists(PATH_now):
		os.makedirs(PATH_now)

	name_energy = PATH_now+str('energy.dat')
	np.savetxt(name_energy, E , fmt='%.9f')

	name_corr = PATH_now+str('corr.dat')
	np.savetxt(name_corr, nn_cor[0] , fmt='%.9f')

quit()















directory = ''
LOCAL = Constants_dictionary.get("LOCAL")

if not os.path.exists(LOCAL+os.sep+directory):
	os.makedirs(LOCAL+os.sep+directory)

PATH_now = LOCAL+os.sep+directory+os.sep

name_fide = PATH_now+str('gap.dat')
np.savetxt(name_fide, lev , fmt='%.9f')

name_fide = PATH_now+str('energy.dat')
np.savetxt(name_fide, E , fmt='%.9f')
















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




