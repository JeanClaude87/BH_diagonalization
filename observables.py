import numpy as np
import os
import time
import function as ff
from scipy.sparse import csc_matrix


## .................................................................
## ....................OBSERVABLES..................................
## .................................................................

#def density():

#..................................................dens
def density(V,**args):
	
	BASE_bose = args.get("BASE_bose")

	den   = np.dot(np.transpose(np.square(np.absolute(V))),BASE_bose)

	return den

def OUTER_creation(A):
	
	Dim = len(A)
	L = len(A[0])
	B = np.zeros((Dim,L,L), dtype=np.float)

	for i in range(Dim):
		B[i] = np.outer(A[i],A[i])
	return B

def NiNj(V,**args):

	states   = args.get("BASE_bose")
	ll  	 = args.get("ll")
	DIM_H 	 = np.int(args.get("DIM_H"))

	Cor_B = np.zeros((DIM_H,ll,ll), dtype=np.float)
	
	for i in range(DIM_H):
		Cor_B[i] = np.outer(states[i],states[i])

	aa = (V.T)[0].T
	NN = np.einsum('n,nij -> ij', np.abs(aa)**2, Cor_B, optimize=True)

	return NN






#############      CdCdCd    CORRELATIONS


def N_creation(**args):

	states   = args.get("BASE_bose")
	ll  	 = np.int(args.get("ll"))
	nn  	 = np.int(args.get("nn"))	
	DIM_H 	 = np.int(args.get("DIM_H"))

	op = np.zeros((ll,DIM_H,DIM_H), dtype=np.float)
	
	for i in range(ll):
			for st in range(DIM_H):	

				op[i,st,st] = states[st,i]

	return op

def NNm1_creation(**args):

	states   = args.get("BASE_bose")
	ll  	 = np.int(args.get("ll"))
	nn  	 = np.int(args.get("nn"))	
	DIM_H 	 = np.int(args.get("DIM_H"))

	op = np.zeros((ll,DIM_H,DIM_H), dtype=np.float)
	
	#print(states)

	for i in range(ll):

		for st in range(DIM_H):	

			op[i,st,st] += states[st,i]*(states[st,i])

	return op

def CdiCj_creation(**args):

	states   = args.get("BASE_bose")
	ll  	 = np.int(args.get("ll"))
	nn  	 = np.int(args.get("nn"))	
	DIM_H 	 = np.int(args.get("DIM_H"))

	mat_type = args.get("mat_type")

	CDC = np.zeros((ll,DIM_H,DIM_H), dtype=np.float)
	
	for i in range(ll-1):
		for st in range(DIM_H):			
			ind, weight = weight_2_ind(i,i+1,st,**args)				
			CDC[i,st,ind] = weight

		print('CDC', i)

	for st in range(DIM_H):	
		ind, weight = weight_2_ind(ll-1,0,st,**args)				
		CDC[ll-1,st,ind] = weight		

	return CDC





def weight_2_ind(i,j,st,**args):

	B_bose 	 = args.get('BASE_bose')	#.......[3 0 0 0 0 0], numpy.ndarray

	peso   = int(1)
	uga    = B_bose[st]*1

	peso   *= uga[j]
	uga[j] -= int(1)

	peso   *= uga[i]+1
	uga[i] += int(1)

	ind = ff.get_index(ff.FROM_bose_TO_bin(uga,**args), **args)	

#	if peso > 0:
	return ind, np.sqrt(peso)
#	else:
#		if peso<0:
#			print(peso)
#		return ind, 0



def CdiCj(V, **args):

	CDC      = args.get("CDC_matrix")
	ll  	 = np.int(args.get("ll"))

	array = np.zeros((ll,ll,4), dtype=np.float)
	
	V_c      = np.conj(V)

	for i in range(ll):	 
		for j in range(ll):

			val = V.dot(CDC[i,j].dot(V_c))

			array[i,j] = [i,j, np.real(val), np.imag(val)] 

	array = array.reshape((ll*ll,4))

	return array



def CdCdCC(V, **args):

	CDC      = args.get("CDC_matrix")
	ll  	 = np.int(args.get("ll"))

	array = np.zeros((ll,ll,ll,ll,6), dtype=np.float)	
	
	V_c	= np.conj(V)

	for i in range(ll):
		for j in range(ll):
			for k in range(ll):
				for l in range(ll):

					val = V.dot(CDC[i,j].dot(CDC[l,k].dot(V_c)))
					array[i,j,l,k] = [i,j,l,k, np.real(val), np.imag(val)]

	array = array.reshape((ll*ll*ll*ll,6))

	return array




def CdiCj_t(psit, Dstep, **args):	

	CDC      = args.get("CDC_matrix")
	ll  	 = np.int(args.get("ll"))
	dt       = args.get("dt")
	step_num = args.get("step_num")
	t_start  = args.get("t_start")

	t_vec 	   = range(0,step_num,Dstep)
	t_num      = len(t_vec)

	array = np.zeros((t_num,ll,ll,5), dtype=np.float)

	for t in t_vec:

		V   = psit[t]
		V_c = np.conj(V)

		for i in range(ll):	 
			for j in range(ll):

				val = V.dot(CDC[i,j].dot(V_c))

				array[t,i,j] = [t*dt+t_start, i, j, np.real(val), np.imag(val)] 

	array = array.reshape((t_num*ll*ll,5))

	return array


def CdCdCC_t(psit, Dstep, **args):

	CDC      = args.get("CDC_matrix")
	ll  	 = np.int(args.get("ll"))
	dt       = args.get("dt")
	step_num = args.get("step_num")
	t_start  = args.get("t_start")

	t_vec 	= range(0,step_num,Dstep)
	t_num 	= len(t_vec)

	array = np.zeros((t_num,ll,ll,ll,ll,7), dtype=np.float)	
	
	for t in t_vec:

		V 	= psit[t]
		V_c	= np.conj(V)

		for i in range(ll):
			for j in range(ll):
				for k in range(ll):
					for l in range(ll):

						val = V.dot(CDC[i,j].dot(CDC[l,k].dot(V_c)))
						array[t,i,j,l,k] = [t*dt+t_start, i, j, l, k, np.real(val), np.imag(val)]

	array = array.reshape((t_num*ll*ll*ll*ll,7))

	return array



def Export_Observable(obs, directory, name, **args):
	
	LOCAL = args.get("LOCAL")

	if not os.path.exists(LOCAL+os.sep+directory):
		os.makedirs(LOCAL+os.sep+directory)

	name_obs = LOCAL+os.sep+directory+os.sep+str(name)
	np.savetxt(name_obs, obs , fmt='%.9f')

	return 0


def Export_Observable_time(psi_t,directory,name,**args):

	dt 	  = args.get("dt")
	nstep = args.get("step_num")
	t_start  = args.get("t_start")
	ll  	 = np.int(args.get("ll"))
	
	DEN   = []

	for i in range(nstep):

		dens = density( psi_t[:,i], **args)

		for j in range(ll):

			DEN.append([i*dt+t_start,j,dens[0,j]])	
	
	Export_Observable(DEN, directory, name, **args)

	return 0


def Export_Fidelity(psi_t, state_B, directory, name,**args):

	dt 	  = args.get("dt")
	nstep = args.get("step_num")
	t_start  = args.get("t_start")

	FID   = []

	for i in range(nstep):

		FID_t = np.square(np.absolute(np.vdot(psi_t[i],state_B)))
		FID.append([i*dt+t_start,FID_t])
	
	Export_Observable(FID, directory, name, **args)

	return 0

def Export_Fidelity_CAT_s(psi_t, psi1, psi2, directory, name,**args):

	ll    = args.get("ll")
	nn    = args.get("nn")
	LOCAL = args.get("LOCAL")
	U  	  = args.get("U")
	bar   = args.get("bar")
	dt 	  = args.get("dt")
	nstep = args.get("step_num")
	t_start  = args.get("t_start")
	
	FID   = []

	for i in range(nstep):

		z1 = np.vdot(psi_t[i],psi1[:,0])
		z2 = np.vdot(psi_t[i],psi2[:,0])
	
		#print(i,np.abs(z1),np.abs(z2))

		zz = (np.conj(z1+z2)*(z1+z2))/2
		
		#print('A', i, zz)

		FID.append([i*dt+t_start, zz, np.conj(z1)*z1 + np.conj(z2)*z2, np.conj(z1)*z2 + np.conj(z2)*z1])

	Export_Observable(np.real(FID), directory, name, **args)

	return 0

def Export_Fidelity_CAT_a(psi_t, psi1, psi2, directory, name,**args):

	ll    = args.get("ll")
	nn    = args.get("nn")
	LOCAL = args.get("LOCAL")
	U  	  = args.get("U")
	bar   = args.get("bar")
	dt 	  = args.get("dt")
	nstep = args.get("step_num")
	t_start  = args.get("t_start")
	
	FID   = []

	for i in range(nstep):

		z1 = np.vdot(psi_t[i],psi1[:,0])
		z2 = np.vdot(psi_t[i],-psi2[:,0])
	
		zz = (np.conj(z1 + z2)*(z1+z2))/2
		
		FID.append([i*dt+t_start,zz])

	Export_Observable(np.real(FID), directory, name, **args)

	return 0


###### HAMILTONIAN OPERATORS

def kinetik_op(omega,**args):

	ll 	   = np.int(args.get("ll"))
	BC 	   = args.get("BC")
	CDC    = args.get("CDC_matrix")
	t      = args.get("t")


	J   = t*np.exp(-2*np.pi*1j*omega/ll)

	op  = (0-0j)*csc_matrix(CDC[0], shape = (DIM_H,DIM_H))

	for i in range(0,ll-1):

		op += csc_matrix(CDC[i], shape = (DIM_H,DIM_H))

	if BC == 0:

		op += csc_matrix(CDC[ll-1], shape = (DIM_H,DIM_H))

	op *= J

	op1 = np.conjugate(op.T)

	return op+op1

def bar_0(x,**args):

	ll 	= np.int(args.get("ll"))
	N = args.get("N_matrix")

	op = csc_matrix(N[x], shape = (DIM_H,DIM_H))
	#op *= float(bar)	

	return op

def corrente_op(omega, **args):

	ll  	 = np.int(args.get("ll"))
	CDC      = args.get("CDC_matrix")
	t 		 = args.get("t") 	
	DIM_H 	 = np.int(args.get("DIM_H"))
	fl  	 = args.get("flux") 

	J   	 = 1j*t*np.exp(-2*np.pi*1j*omega/ll)

	BC 	   	 = args.get("BC")


	op  = (0-0j)*csc_matrix(CDC[0], shape = (DIM_H,DIM_H))

	for i in range(0,ll-1):

		op += csc_matrix(CDC[i], shape = (DIM_H,DIM_H))

	if BC == 0:

		op += csc_matrix(CDC[ll-1], shape = (DIM_H,DIM_H))

	op *= J

	op1 = np.conjugate(op.T)

	return op+op1

def fluct_op(op,**args):

	fl_op = np.matmul(op,op)
    
	return fl_op



def int_op(**args):
	
	ll  	 = np.int(args.get("ll"))
	N        = args.get("N_matrix")
	NN       = args.get("NN_matrix")

	Hint = 0*csc_matrix(N[0], shape = (DIM_H,DIM_H))
	
	for j in range(ll):
		Hint += csc_matrix(NN[j] - N[j], shape = (DIM_H,DIM_H))

	return Hint




























