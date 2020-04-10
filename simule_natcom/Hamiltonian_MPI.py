import numpy as np
from scipy.sparse import csc_matrix
from scipy.sparse import linalg as linalgS
from numpy import linalg as lin
import time

import hamiltonian        as ham
import hamiltonian_parity as ham_par
import function           as ff
import observables        as ob
import time_evolution	  as t_ev

from mpi4py import MPI

def split(container, count):
	"""
	Simple function splitting a container into equal length chunks.
	Order is not preserved but this is potentially an advantage depending on
	the use case.
	"""
	return [container[_i::count] for _i in range(count)]


def bose_Hamiltonian ( **args):


#Hamiltonian needs as input:


#...... Parameter: DIM_H
	if COMM.rank == 0:

		DIM_H 	 = np.int(args.get("DIM_H"))
		ll 		 = args.get("ll")	
		nn 		 = args.get("nn")		

	#Hamiltonian returns:
	#...... Sparse or Dense matrix. By default is SPARSE
		
		mat_type = args.get("mat_type")
		if mat_type == None:
			mat_type = 'Sparse'

	
	############### MPI VERSION

	if COMM.rank == 0:
		jobs = list(range(DIM_H))
		jobs = split(jobs, COMM.size)
	else:
		jobs = None

	jobs = COMM.scatter(jobs, root=0)

	XX = []
	YY = []
	AA = []

	for i in jobs:
		res = ham.evaluate_ham(i, **args)
		XX.append(res[0])
		YY.append(res[1])
		AA.append(res[2])

	XX0 = MPI.COMM_WORLD.gather( XX, root=0)
	YY0 = MPI.COMM_WORLD.gather( YY, root=0)
	AA0 = MPI.COMM_WORLD.gather( AA, root=0)

	if COMM.rank == 0:

		X0 = [item for sublist in XX0 for item in sublist]
		Y0 = [item for sublist in YY0 for item in sublist]
		A0 = [item for sublist in AA0 for item in sublist]

		print("Results:", 'porcodio')

		X1 = [item for sublist in X0 for item in sublist]
		Y1 = [item for sublist in Y0 for item in sublist]
		A1 = [item for sublist in A0 for item in sublist]


		Hamiltonian = csc_matrix((A1, (X1,Y1)), shape=(DIM_H,DIM_H), dtype=np.double)
		ff.print_matrix(Hamiltonian)

		if mat_type == 'Dense':

			Hamiltonian = csc_matrix.todense(Hamiltonian)



		return Hamiltonian
