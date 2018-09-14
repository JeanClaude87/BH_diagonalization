import numpy as np
from math import factorial
import math
import itertools
from scipy.sparse import csc_matrix
from scipy.sparse import linalg as linalgS
from numpy import linalg as lin
from numpy import matlib








#################################################################
####   For comparing different approaches to diagonalize     ####
#################################################################



#~ np.set_printoptions(precision=2)
#~ PATH_now = os.path.abspath('.')


#~ #........ LIST OF GLOBAL PARAMETERS

#~ ### -> ll, nn, tab_fact, DIM_H, BASE_bin, BASE_bose, CORR_BASE

#~ ll           = 8
#~ ff.ll        = ll

#~ nn           = 6
#~ ff.nn        = nn

#~ ff.tab_fact  = tab_fact   = ff.fact_creation(nn+ll)

#~ ff.DIM_H     = DIM_H      = ff.hilb_dim(nn,ll)

#~ BASE_bin, BASE_bose       = ff.Base_prep()

#~ ff.BASE_bin  = BASE_bin
#~ ff.BASE_bose = BASE_bose

#~ CORR_BASE    = ff.OUTER_creation(BASE_bose)
#~ ff.CORR_BASE = CORR_BASE


#~ #........ LIST OF HAMILTONIAN PARAMETERS

#~ t=-1.
#~ U=-3.0

#~ BC=0




#~ nstate = DIM_H

#~ base_parity_ind = ff.base_parity()

#print(base_parity_ind)

#~ DIM_par_H = len(base_parity_ind)



#~ ################################################ Original Hamiltonian

#~ print('Preparing Base')
#~ time0=time.clock()
#~ ham_ind1, ham_ind2, ham_val = ff.bose_Hamiltonian(BC,t,U)
#~ time1=time.clock()
#~ print ('')
#~ print ('Time=',time0,time1,"%e" % (time1-time0))
#~ print ('')


#print('original bose base',BASE_bose)
#~ Sp_Hamiltonian = ff.make_sparse_mat(ham_ind1, ham_ind2, ham_val, DIM_H)

#ff.print_hamiltonian(Sp_Hamiltonian)

#~ print('Diagonalizing')
#~ time0=time.clock()
#~ ED1 ,EV1 = ff.diagonalization(Sp_Hamiltonian,DIM_H-1,True)
#~ time1=time.clock()
#~ print ('')
#~ print ('Time=',time0,time1,"%e" % (time1-time0))
#~ print ('')

#print(ED1)


#~ ################################################ PIERO STUFF

#~ print('Preparing Base')
#~ time0=time.clock()
#~ H_par = ff.bose_Hamiltonian_parity(Sp_Hamiltonian,base_parity_ind,BC,t,U,1)
#~ time1=time.clock()
#~ print ('')
#~ print ('Time=',time0,time1,"%e" % (time1-time0))
#~ print ('')

#ff.print_hamiltonian(H_par)

#~ print('Diagonalizing')
#~ time0=time.clock()
#~ ED ,EV = ff.diagonalization(H_par,DIM_H-1,True)
#~ time1=time.clock()
#~ print ('')
#~ print ('Time=',time0,time1,"%e" % (time1-time0))
#~ print ('')


#print(ED)



#~ ################################################ My stuff

#~ print('Preparing Base')
#~ time0=time.clock()
#~ ham_p_ind1, ham_p_ind2, ham_p_val = ff.bose_Hamiltonian_parity2(base_parity_ind,ham_ind1,ham_ind2,ham_val,BC,t,U,1)
#~ time1=time.clock()
#~ print ('')
#~ print ('Time=',time0,time1,"%e" % (time1-time0))
#~ print ('')


#~ SP_Hamiltonian = ff.make_sparse_mat(ham_p_ind1, ham_p_ind2, ham_p_val, DIM_H)

#ff.print_hamiltonian(SP_Hamiltonian)

#~ print('Diagonalizing')
#~ time0=time.clock()
#~ ED ,EV = ff.diagonalization(SP_Hamiltonian,2,False)
#~ time1=time.clock()
#~ print ('')
#~ print ('Time=',time0,time1,"%e" % (time1-time0))
#~ print ('')

#print(ED)

























def bose_Hamiltonian_parity2(base_parity_ind,ham_ind1,ham_ind2,ham_val,BC,t,U,parity):

	p=parity

	ham_p_ind1 = []
	ham_p_ind2 = []
	ham_p_val  = []


	for i in range(len(base_parity_ind)):

		x = base_parity_ind[i]


##----- INTERACTIONS

		if x[0] == x[1]:
			if p == -1:
				continue
			else:
				for j in range(len(ham_ind1)):
					if x[0] == ham_ind1[j]:												### This is the 0th element of the state of the parity basis
						#~ print(x[0],ham_ind2[j],ham_val[j])
	
						for k in range(len(base_parity_ind)):
							y = base_parity_ind[k]
							if y[0] == ham_ind2[j] or y[1] == ham_ind2[j]:
								ham_p_ind1.append( i )
								ham_p_ind2.append( k )
								if y[0]==y[1]:
									ham_p_val.append( ham_val[j] )	
								else:
									ham_p_val.append( (1/np.sqrt(2))*ham_val[j] )
								#~ print('Final Indeces','I=',i,'J=',k,'value=',ham_p_val[k])
										

		else:
			hola=0
			for j in range(len(ham_ind1)):												### We run over all possible indices in our sparse Hamiltonian
				for g in range(2):
					if x[g] == ham_ind1[j]:												### We need to look for all connection of the composite state
						#~ print(x[0],ham_ind2[j],ham_val[j])
	
						for k in range(len(base_parity_ind)):
							y = base_parity_ind[k]
							if y[0] == ham_ind2[j] or y[1] == ham_ind2[j]:				### Goes through all the possible couplings existing between state x[0] and the rest
								ham_p_ind1.append( i )
								ham_p_ind2.append( k )
								if y[0]==y[1]:   										### To check that we are dealing with 1 state formed by 2
									ham_p_val.append( (1/np.sqrt(2))*ham_val[j] )	
								else:													### To check that we are dealing with 2 states formed by 2
									ham_p_val.append( (1/2)*ham_val[j] )				
								#~ print('Final Indeces','I=',i,'J=',k,'value=',ham_p_val[k])
													
									
				
#~ # Parity -1
	index_p_neg=len(base_parity_ind)
	index_p_neg2=len(base_parity_ind)
	for i in range(len(base_parity_ind)):
		
		x = base_parity_ind[i]														#### This selects the first parity state
		if x[0] == x[1]:															#### We only care about composite states
			index_p_neg -=1
			continue
		else:
			for j in range(len(ham_ind1)):	
				for g in range(2):													#### We go through the two elements of the composite state
					if x[g] == ham_ind1[j]:											#### We take the elements of the original Hamiltonian that start from the first element of the composite state					
						#~ print('')
						#~ print(x[g],ham_ind2[j],ham_val[j])						
						index_p_neg2=len(base_parity_ind)							#### To avoid jumping of indices when going through a non composite state
						for k in range(len(base_parity_ind)):						
							y = base_parity_ind[k]									#### We run over all the states of the parity basis
							if y[0]!=y[1]:											#### This checks that we only consider composed states
								if y[0] == ham_ind2[j] or y[1] == ham_ind2[j]:		#### If there exist the possibility to jump to other states we take those possibilities
									if x[g]!=ham_ind2[j]:
										if g==0 and y[1]==ham_ind2[j] or g==1 and y[0]==ham_ind2[j]:
											value=-(1/2)*ham_val[j]
											ham_p_ind1.append( i+index_p_neg )
											ham_p_ind2.append( k+index_p_neg2 )
											ham_p_val.append( value )									
											#~ print('g=',g,'x[g]=',x[g],'index[j]=',ham_ind2[j],value)										
										else:
											value=(1/2)*ham_val[j]
											ham_p_ind1.append( i+index_p_neg )
											ham_p_ind2.append( k+index_p_neg2 )
											ham_p_val.append( value )
											#~ print('g=',g,'x[g]=',x[g],'index[j]=',ham_ind2[j],value)																													
									else:	
										value=(1/2)*ham_val[j]
										ham_p_ind1.append( i+index_p_neg )
										ham_p_ind2.append( k+index_p_neg2 )
										ham_p_val.append( value )	
										#~ print('g=',g,'x[g]=',x[g],'index[j]=',ham_ind2[j],value)									
								
								
																
							else:
								index_p_neg2 -=1
				

	return  ham_p_ind1,ham_p_ind2,ham_p_val

