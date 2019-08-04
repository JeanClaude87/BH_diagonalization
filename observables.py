import numpy as np
import os
import time
import function as ff
from scipy import sparse

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

def CdiCj_creation(**args):

	states   = args.get("BASE_bose")
	ll  	 = np.int(args.get("ll"))
	nn  	 = np.int(args.get("nn"))	
	DIM_H 	 = np.int(args.get("DIM_H"))

	B_bose 	 = args.get('BASE_bose')	#.......[3 0 0 0 0 0], numpy.ndarray

	Cd  = np.power(np.remainder(states+1,nn+1),1/2)
	C   = np.power(states,1/2)

	CDC = np.zeros((ll,ll,DIM_H,DIM_H), dtype=np.float)
	
	for i in range(ll):
		for j in range(ll):

			if i!=j:
				for st in range(DIM_H):		

					if(Cd[st,i]*C[st,j] > 0):

						uga    = B_bose[st]*1
						uga[i] += 1
						uga[j] -= 1

						#state_hop = B_bose[st]+hops
						ind = ff.get_index(ff.FROM_bose_TO_bin(uga,**args), **args)	

						weight = np.sqrt((B_bose[st,i]+1)*B_bose[st,j])
								 #np.sqrt((B_bose[st,i]+1)*B_bose[st,j])
						
						CDC[i,j,st,ind] = weight

						#print('robina', 'create', i+1, 'destroy', j+1, B_bose[st], uga, weight)
				#if B_bose[st,j] <= 0: print('non puoi distruggere')
				#if B_bose[st,i] >= 4: print('non puoi creare')

	return CDC


def CdCdCC_creation(**args):

	ll  	 = np.int(args.get("ll"))
	nn  	 = np.int(args.get("nn"))	
	DIM_H 	 = np.int(args.get("DIM_H"))

	B_bose 	 = args.get('BASE_bose')	#.......[3 0 0 0 0 0], numpy.ndarray

	CDC = np.zeros((ll,ll,ll,ll,DIM_H,DIM_H), dtype=np.float)
	Cd  = np.power(np.remainder(B_bose+1,nn+1),1/2)
	C   = np.power(B_bose,1/2)	
	
	for i in range(ll):
		for j in range(ll):
			for k in range(ll):
				for l in range(ll):
					for st in range(DIM_H):	

						ind, weight = weight_4_ind(i,j,k,l,st,**args)
						CDC[i,j,k,l,st,ind] = weight	

	return CDC

def weight_4_ind(i,j,k,l,st,**args):

	ll  	 = np.int(args.get("ll"))
	nn  	 = np.int(args.get("nn"))	
	B_bose 	 = args.get('BASE_bose')	#.......[3 0 0 0 0 0], numpy.ndarray

	peso   = 1

	uga    = B_bose[st]*1

	peso   *= uga[l]
	uga[l] -= 1

	peso   *= uga[k]+1
	uga[k] += 1

	peso   *= uga[j]
	uga[j] -= 1

	peso   *= uga[i]+1
	uga[i] += 1

	ind = ff.get_index(ff.FROM_bose_TO_bin(uga,**args), **args)	

	if peso > 0:
		return ind, np.sqrt(peso)
	else:
		return ind, 0

def CdCdCC(V, **args):

	CDCDCC 	= args.get("CDCDCC_matrix")
	V_c  	= np.conj(V)

	correl  = np.einsum('xyztlj,l,j -> xyzt', CDCDCC, V, V_c) #, optimize=True)

	return correl

def CdiCj(V, **args):

	CDC      = args.get("CDC_matrix")
	
	V_c      = np.conj(V)
	dens 	 = density(V, **args)

	CdiCj  = np.einsum('xylj,l,j -> xy', CDC, V, V_c, optimize=True)
	CdiCj += np.diag(dens)

	return CdiCj

def CdCdCC_t(psit, Dstep, **args):

	CDCDCC 	= args.get("CDCDCC_matrix")

	ll  	 = np.int(args.get("ll"))
	dt       = args.get("dt")
	step_num = args.get("step_num")

	t_vec 	   = range(0,step_num,Dstep)
	prop_array = np.array([[[[[[t*dt, i, j, k, l, psit[t].dot(CDCDCC[i,j,k,l].dot(np.conj(psit[t])))] for t in t_vec] for i in range(ll)] for j in range(ll)] for k in range(ll)] for l in range(ll)]).reshape((step_num*ll*ll*ll*ll,6))

	print(prop_array.shape)

	return prop_array


def CdiCj_t(psit, Dstep, **args):	


	CDC      = args.get("CDC_matrix")

	ll  	 = np.int(args.get("ll"))
	dt       = args.get("dt")
	step_num = args.get("step_num")

	t_vec 	   = range(0,step_num,Dstep)
	prop_array = np.array([[[[t*dt, i, j, psit[t].dot(CDC[i,j].dot(np.conj(psit[t])))] for t in t_vec] for i in range(ll)] for j in range(ll)] ).reshape((step_num*ll*ll,4))

	return prop_array


def corrente(V, **args):

	ll  	 = np.int(args.get("ll"))

	CdC 	 = CdiCj(V, **args)
	xx 		 = -np.imag([(CdC[i,i+1]-CdC[i+1,i]) for i in range(ll-1) ])
	corrente = np.sum(xx)

	return corrente

def corrente_t(psit, Dstep, **args):

	ll  	 = np.int(args.get("ll"))
	dt       = args.get("dt")
	step_num = args.get("step_num")
	
	t_vec          = range(0,step_num,Dstep)
	corrente_array = np.array([[t*dt, corrente(psit[t], **args)] for t in t_vec])

	return corrente_array


def Olsh2(V, **args):

	states   = args.get("BASE_bose")
	ll  	 = args.get("ll")
	nn  	 = args.get("nn")
	DIM_H 	 = np.int(args.get("DIM_H"))

	Cor_B = np.zeros((DIM_H,ll,ll), dtype=np.float)
	
	coeff 	= [[ (i-j)**2 for i in range(ll)] for j in range(ll)]
	coeff	= np.asmatrix(coeff)

	for i in range(DIM_H):
		Cor_B[i] = np.outer(states[i],states[i])

	aa = (V.T)[0].T
	ol2 = np.einsum('n,nij,ij -> ij', np.abs(aa)**2, Cor_B, coeff, optimize=True)/(nn**2)

	return ol2

def Olsh1(V, **args):

	states   = args.get("BASE_bose")
	ll  	 = args.get("ll")
	nn  	 = args.get("nn")
	DIM_H 	 = np.int(args.get("DIM_H"))

	BASE_bose = args.get("BASE_bose")

	den   = np.dot(np.transpose(np.square(np.absolute(V))),BASE_bose)

	coeff 	= [ den[1,i]*(i-3*ll/4)**2 for i in range(ll)]

	ol1 = np.sum(coeff)/(nn)

	return ol1


def Export_Observable(obs, directory, name, **args):
	
	LOCAL = args.get("LOCAL")

	if not os.path.exists(LOCAL+os.sep+directory):
		os.makedirs(LOCAL+os.sep+directory)

	name_obs = LOCAL+os.sep+directory+os.sep+str(name)
	np.savetxt(name_obs, obs , fmt='%.9f')

	return 0


def Export_Observable_time(psi_t,directory,name,**args):

	ll    = args.get("ll")
	nn    = args.get("nn")
	LOCAL = args.get("LOCAL")
	U  	  = args.get("U")
	dt 	  = args.get("dt")
	nstep = args.get("step_num")
	
	DEN   = []

	for i in range(nstep):

		dens = density( psi_t[:,i], **args)

		for j in range(ll):

			DEN.append([i*dt,j,dens[0,j]])	
	
	Export_Observable(FID, directory, name, **args)

	return 0


def Export_Fidelity(psi_t, state_B, directory, name,**args):

	ll    = args.get("ll")
	nn    = args.get("nn")
	LOCAL = args.get("LOCAL")
	U  	  = args.get("U")
	bar   = args.get("bar")
	dt 	  = args.get("dt")
	nstep = args.get("step_num")

	FID   = []

	for i in range(nstep):

		FID_t = np.square(np.absolute(np.vdot(psi_t[i],state_B)))
		FID.append([i*dt,FID_t])
	
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
	
	FID   = []

	for i in range(nstep):

		z1 = np.vdot(psi_t[i],psi1[:,0])
		z2 = np.vdot(psi_t[i],psi2[:,0])
	
		zz = np.real(( z1*np.conj(z1) + z2*np.conj(z2) + z1*np.conj(z2) + z2*np.conj(z1) ) / 2.0 )
		FID.append([i*dt,zz])

	Export_Observable(FID, directory, name, **args)

	return 0

def Export_Fidelity_CAT_a(psi_t, psi1, psi2, directory, name,**args):

	ll    = args.get("ll")
	nn    = args.get("nn")
	LOCAL = args.get("LOCAL")
	U  	  = args.get("U")
	bar   = args.get("bar")
	dt 	  = args.get("dt")
	nstep = args.get("step_num")
	
	FID   = []

	for i in range(nstep):

		z1 = np.vdot(psi_t[i],psi1[:,0])
		z2 = np.vdot(psi_t[i],psi2[:,0])
	
		zz = np.real(( z1*np.conj(z1) + z2*np.conj(z2) - z1*np.conj(z2) - z2*np.conj(z1) ) / 2.0 )
		FID.append([i*dt,zz])

	Export_Observable(FID, directory, name, **args)

	return 0





