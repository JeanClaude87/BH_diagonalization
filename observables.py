import numpy as np
import os

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


def CdiCj(V, dens, **args):

	states   = args.get("BASE_bose")
	ll  	 = args.get("ll")
	nn  	 = args.get("nn")	

	Cd  = np.power(np.remainder(states+1,nn+1),1/4)
	C   = np.power(states,1/4)
	V_c = np.conj(V)

	CdiCj = np.zeros((ll,ll), dtype=np.complex)

	for i in range(ll):
		for j in range(i+1,ll):

				Cd_C = Cd[:,i] * C [:,j]	# V
				C_Cd = C [:,i] * Cd[:,j] 	# V*

				A = V   * Cd_C
				B = V_c * C_Cd
				
				CdiCj[i,j] = np.dot(A[A!=0],B[B!=0])

				
	CdiCj += CdiCj.T
	CdiCj += np.diag(dens)

	return CdiCj


def Export_Observable():

	a = 0

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


def Export_Fidelity(psi_t,dt,name,**args):

	ll    = args.get("ll")
	nn    = args.get("nn")
	LOCAL = args.get("LOCAL")
	U  	  = args.get("U")
	
	nstep = len(psi_t.T)
	FID   = []

	A = np.squeeze(np.asarray(psi_t[:,0]))

	for i in range(nstep):

		B = np.squeeze(np.asarray(psi_t[:,i]))
		FID.append([i*dt,np.square(np.absolute(np.dot(A,B)))])

	directory = '/dati/L_'+str(ll)+'-N_'+str(nn)+os.sep+'U_'+str(U)
	
	if not os.path.exists(LOCAL+os.sep+directory):
		os.makedirs(LOCAL+os.sep+directory)

	PATH_now = LOCAL+os.sep+directory+os.sep

	name_fide = PATH_now+str(name)
	np.savetxt(name_fide, FID , fmt='%.9f')

	return 0







