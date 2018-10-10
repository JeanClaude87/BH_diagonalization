import numpy as np

import hamiltonian        as ham
import hamiltonian_parity as ham_par
import function           as ff
import observables        as ob



def inital_state(part_ind,**args):

	b_p_inp	 = args.get("parity_index")
	b_p      = np.asarray(b_p_inp)	
	DX       = args.get("sim_sec_len")

	PARITY   = args.get("parity")

	nn 		 = args.get("nn")
	ll 		 = args.get("ll")
	BASE_bose = args.get("BASE_bose")
	DIM_H    = args.get("DIM_H")

	state    = np.zeros(ll, dtype=np.int)

	for x in range(nn):
		state[part_ind[x]] += int(1)

	state_con = ff.FROM_bose_TO_bin (state,     **args)
	state_ind = ff.get_index        (state_con, **args)
	B 		  = np.zeros(DIM_H, dtype=np.double)
	
	if  PARITY == 'True':

		state_rev_ind  = ham_par.parity      (state_con, **args)[1]

		ind      = min(state_ind,state_rev_ind)
		par_ind  = b_p[ind]
		
		i_s = par_ind[0]
		i_a = par_ind[3]+DX

		if   state_ind == state_rev_ind:

			B[i_s] = 1

		elif state_ind < state_rev_ind:

			B[i_s] = +np.sqrt(2)/2
			B[i_a] = +np.sqrt(2)/2

		elif state_ind > state_rev_ind:

			B[i_s] = +np.sqrt(2)/2
			B[i_a] = -np.sqrt(2)/2

	else:

		B[state_ind] += 1	 

	return B



