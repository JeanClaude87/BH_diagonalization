import numpy as np
import os
import function as ff

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
	NN = np.einsum('n,nij -> ij', np.abs(aa)**2, Cor_B)

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
		for j in range(i+1,ll):

			hop = np.zeros(ll, dtype=np.int)

			hop[i] = int(1)
			hop[j] = -1*int(1)

			for st in range(DIM_H):		

				if(Cd[st,i]*C[st,j] > 0):

					state_hop = B_bose[st]+hop
					ind = ff.get_index(ff.FROM_bose_TO_bin(state_hop,**args), **args)	

					weight = np.sqrt((B_bose[st,i]+1)*B_bose[st,j])
					
					CDC[i,j,st,ind] = weight

			CDC[i,j] += CDC[i,j].T

	return CDC

def CdiCj(V, **args):

	CDC   = args.get("CDC_matrix")
	flux 	 = np.int(args.get("U"))
	ll  	 = np.int(args.get("ll"))
	
	V_c = np.conj(V)
	dens 	 = density(V, **args)

	#CdiCj

	CdiCj = np.einsum('xylj,l,j -> xy', CDC, V_c, V)

	CdiCj += CdiCj.T
	CdiCj *= 0.5
	CdiCj += np.diag(dens)

	for ii in range(ll):
		for jj in range(ll):
			CdiCj *= np.exp(2*np.pi*1j*ii*flux/ll)*np.exp(-2*np.pi*1j*jj*flux/ll)


	return CdiCj




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
	ol2 = np.einsum('n,nij,ij -> ij', np.abs(aa)**2, Cor_B, coeff)/(nn**2)

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


def Export_Observable_time(psi_t,dt,name,**args):

	ll    = args.get("ll")
	nn    = args.get("nn")
	LOCAL = args.get("LOCAL")
	U  	  = args.get("U")
	
	nstep = len(psi_t.T)
	DEN   = []

	for i in range(nstep):

		dens = density( psi_t[:,i], **args)

		for j in range(ll):

			DEN.append([i*dt,j,dens[0,j]])	

	directory = '/dati/L_'+str(ll)+'-N_'+str(nn)+os.sep+'U_'+str(U)
	
	if not os.path.exists(LOCAL+os.sep+directory):
		os.makedirs(LOCAL+os.sep+directory)

	PATH_now = LOCAL+os.sep+directory+os.sep

	name_dens = PATH_now+str(name)
	np.savetxt(name_dens, DEN , fmt='%.9f')

	return 0


def Export_Fidelity(psi_t,state_B,dt,name,**args):

	ll    = args.get("ll")
	nn    = args.get("nn")
	LOCAL = args.get("LOCAL")
	U  	  = args.get("U")
	bar   = args.get("bar")
	
	nstep = len(psi_t)

	FID   = []

	for i in range(nstep):
		FID_t = np.square(np.absolute(np.vdot(psi_t[i],state_B)))
		FID.append([i*dt,FID_t])

	directory = os.sep+'dati'+os.sep+'L_'+str(ll)+os.sep+'N_'+str(nn)+os.sep+'U_'+str(U)+os.sep+'bb_'+str(bar)
#	directory = '/dati/L_'+str(ll)+'/N_'+str(nn)+os.sep+'U_'+str(U)
	
	if not os.path.exists(LOCAL+os.sep+directory):
		os.makedirs(LOCAL+os.sep+directory)

	PATH_now = LOCAL+os.sep+directory+os.sep

	name_fide = PATH_now+str(name)
	np.savetxt(name_fide, FID , fmt='%.9f')

	'''
	FID   = []

	A = np.squeeze(np.asarray(psi_t[:,0]))


		B = np.squeeze(np.asarray(psi_t[:,i]))
		FID.append([i*dt,np.square(np.absolute(np.dot(A,B)))])

		print(i*dt)

	directory = '/dati/L_'+str(ll)+'/N_'+str(nn)+os.sep+'U_'+str(U)
	
	if not os.path.exists(LOCAL+os.sep+directory):
		os.makedirs(LOCAL+os.sep+directory)

	PATH_now = LOCAL+os.sep+directory+os.sep

	name_fide = PATH_now+str(name)
	np.savetxt(name_fide, FID , fmt='%.9f')
	'''

	return 0

def Export_Fidelity_time(psi_t,dt,name,**args):

	ll    = args.get("ll")
	nn    = args.get("nn")
	LOCAL = args.get("LOCAL")
	U  	  = args.get("U")
	bar   = args.get("bar")
	
	nstep = len(psi_t.T)
	FID   = []

	A = np.squeeze(np.asarray(psi_t[:,0]))

	for i in range(nstep):

		B = np.squeeze(np.asarray(psi_t[:,i]))
		FID.append([i*dt,np.square(np.absolute(np.dot(A,B)))])

		#print(i*dt)

	directory = os.sep+'dati'+os.sep+'L_'+str(ll)+os.sep+'N_'+str(nn)+os.sep+'U_'+str(U)+os.sep+'bb_'+str(bar)
	print(directory)

	if not os.path.exists(LOCAL+os.sep+directory):
		os.makedirs(LOCAL+os.sep+directory)

	PATH_now = LOCAL+os.sep+directory+os.sep

	name_fide = PATH_now+str(name)
	np.savetxt(name_fide, FID , fmt='%.9f')

	return 0







