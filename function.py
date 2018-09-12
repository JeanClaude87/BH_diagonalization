import numpy as np
from math import factorial
import math
import itertools
from scipy.sparse import csc_matrix
from scipy.sparse import linalg as linalgS
from numpy import linalg as linalgD

#..................................hilbert space dimension
def fact_creation(a):
	tab = [factorial(x) for x in range(a)]

	#...... put the faster way

	return tab

def hilb_dim(x,y):

	kk = tab_fact[x+y-1]/(tab_fact[x]*tab_fact[y-1])
	uga= int(kk)
	
	return uga

#..................................base preparation
def Base_prep():
		base_bin = []
		base_num = []
		for bits in itertools.combinations(range(nn+ll-1), nn):
			s = ['0'] * (nn+ll-1)	
			for bit in bits:
				s[bit] = '1'		
			base_bin.append(''.join(s))

			bose = TO_bose_conf(s)
			base_num.append(bose)

		base_bose = np.asarray(base_num, dtype=np.int8)

		return base_bin, base_bose


#..................................parity transformation
# A states
def parity(state):
	
	parity_state 	= state[::-1]

	index_ps = get_index(parity_state)
	
	return parity_state,index_ps

def base_parity():

	base_par = []
	UGA = [i for i in range(DIM_H)]

	for i in range(DIM_H):

		if UGA[i] == "hola" :
			continue

		p_state,index_p_state = parity(BASE_bin[i])

		UGA[index_p_state] 	  = "hola"
		
		base_par.append((i,index_p_state))

	return base_par


#..................................index search
def get_index(state):

	size   = int(nn+ll-1)
	r_par  = int(nn)    #remaining_particles
	r_sit  = int(ll)	#remaining_sites
	result = DIM_H

	#print(size, len(state), state)

	for jj in range(size):
		action_i = int(state[jj]);
		if r_par==0:
			break
		if action_i == 0:
			#print(jj,r_par)		
		
			result -= hilb_dim(r_par-1, r_sit);
			r_sit  -= 1;

		else:
			r_par -= 1;
			#print('else',jj,r_par)

	return DIM_H-result

#..................................from configuration to bin number
def TO_bin(xx):
	return int(xx,2)

#..................................from bin number to configuration
def TO_con(x,L):
	x1=int(x)
	L1=int(L)
	return np.binary_repr(x1, width=L1)

#..................................hop. preparation
#............. ......BC=0 periodic, BC=1 open
def Hop_prep(X,Y):

	Hop_dim=X+Y-2
	
	return [TO_con(2**i+2**((i+1)%(X+Y-1)),X+Y-1) for i in range(Hop_dim)]



#..................................counting number of one
POPCOUNT_TABLE16 = [0] * 2**16
for index in range(len(POPCOUNT_TABLE16)):
	POPCOUNT_TABLE16[index] = (index & 1) + POPCOUNT_TABLE16[index >> 1]

def one_count(v):
	return (POPCOUNT_TABLE16[ v        & 0xffff] +
			POPCOUNT_TABLE16[(v >> 16) & 0xffff])



#from 00100111001 to 010301 (or come cazzo si fa)
def TO_bose_conf(x):
	p=0
	conf = np.zeros(ll, dtype=np.int)
	for jj in range(len(x)):
		if x[jj]=='1':
			conf[p]+=1
		else: 
			p+=1	
	return 	conf


def kin_expval(x,y):

	state_x   = BASE_bin[x]
	bosecon_x = TO_bose_conf(state_x)

	state_y   = BASE_bin[y]
	bosecon_y = TO_bose_conf(state_y)

	jump_c  =	np.argmax(bosecon_x-bosecon_y)
	jump_cd  =	np.argmin(bosecon_x-bosecon_y)

	result	  = np.sqrt((bosecon_x[jump_cd]+1)*(bosecon_x[jump_c])) 
	
	return result

def action_interactions(state,U):
	bosecon = TO_bose_conf(state)

	int_val = 0
	for x in range(ll):
		nx = bosecon[x]
		int_val += U*nx*(nx-1.)/(2.) #+10**-3*np.random.random()

	return int_val 

def bose_Hamiltonian (BC,t,U):

	ham_ind1 = []
	ham_ind2 = []
	ham_val  = []

	for i in range(DIM_H):
		state = BASE_bin[i]

##----- INTERACTIONS

		int_val = action_interactions(state,U)

		#---- INTERACTION = we store i,i,int_val !!!!
		ham_ind1.append( i )
		ham_ind2.append( i )
		ham_val.append( int_val )



##----- KINETIC
		for hop in Hop_prep(ll,nn):
			hop_state_bin = TO_bin(state)^TO_bin(hop)
			hop_state     = TO_con(hop_state_bin,ll+nn-1)
			
			# we cut states with not N particles
			if one_count(hop_state_bin) == nn:

				j = get_index(hop_state)	
				kin_val = t*kin_expval(i,j)

		#---- KINETIC = we store i,j,kin_val !!!!
				ham_ind1.append( i )
				ham_ind2.append( j )
				ham_val.append( kin_val )
	
		# here we put the PERIODIC
		#...................BC=0 periodic, BC=1 open

		if BC == 0:
			if state[0] == '1':

				PBC_newstate1 = state[1:]+state[0]				
				j = get_index(PBC_newstate1)	

				kin_val = t*kin_expval(i,j)
				
		#---- KINETIC = we store i,j,kin_val !!!!
				ham_ind1.append( i )
				ham_ind2.append( j )
				ham_val.append( kin_val )

			if state[-1] == '1':

				PBC_newstate2 = state[-1]+state[:-1]
				j = get_index(PBC_newstate2)	

				kin_val = t*kin_expval(i,j)
				
		#---- KINETIC = we store i,j,kin_val !!!!
				ham_ind1.append( i )
				ham_ind2.append( j )
				ham_val.append( kin_val )
		
	return ham_ind1,ham_ind2,ham_val


# 

def bose_Hamiltonian_parity(base_parity_ind,BC,t,U,parity):

	p=parity

	ham_ind1 = []
	ham_ind2 = []
	ham_val  = []

	for i in range(len(base_parity_ind)):

		x = base_parity_ind[i]

		print(x)

##----- INTERACTIONS

		if x[0] == x[1]:
			if p == -1:
				continue
			else:
				int_val = action_interactions(BASE_bin[x[0]],U)

		else:		
			int_val =   (1/np.sqrt(2))*(action_interactions(BASE_bin[x[0]],U))
			int_val +=	(p/np.sqrt(2))*(action_interactions(BASE_bin[x[1]],U))

		#---- INTERACTION = we store i,i,int_val !!!!
#		ham_ind1.append( i )
#		ham_ind2.append( i )
#		ham_val.append( int_val )

	return  ham_ind1,ham_ind2,ham_val


def diagonalization(X,Y,A_XY,DIM,num_eig,sparse):

# X 		-> vector index i
# Y 		-> vector index j
# A_XY 		-> matrix element A[i,j]
# num_eig 	-> how many eigenvalues: less than DIM_H

	if num_eig >= DIM_H:
		num_eig = DIM_H-4

	numpy_ind1 = np.asarray(X)
	numpy_ind2 = np.asarray(Y)
	numpy_val  = np.asarray(A_XY)


	Hamiltonian = csc_matrix((numpy_val, (numpy_ind1, numpy_ind2)), shape=(DIM,DIM), dtype=np.double)
	#...... tol=10**-20

	if sparse == True:

		eig = linalgS.eigsh(Hamiltonian, k=num_eig, which='SA', return_eigenvectors=True)

	else:

		hamdens = csc_matrix.todense(Hamiltonian)
		eig  = linalgD.eigh(hamdens)
		eigl = list(eig)
		eigl[1] = np.squeeze(np.asarray(eig[1]))
		eig = eigl

	return eig




## .................................................................
## ....................OBSERVABLES..................................
## .................................................................

#def density():

#..................................................dens
def density(V):

	den   = np.dot(np.transpose(V**2),Base_bose)

	return den

def OUTER_creation(A):
	
	Dim = len(A)
	L = len(A[0])
	B = np.zeros((Dim,L,L), dtype=np.float)

	for i in range(Dim):
		B[i] = np.outer(A[i],A[i])
	return B

def NiNj(V,matrix):

	NN = np.einsum('n,nij -> ij', V**2, matrix)

	return NN

def NfixNr(i,V,matrix):

	NiN = np.einsum('n,nj -> j', V**2, matrix[:,i])

	return NiN



#..................................................generate filename
def generate_filename(basename):
	unix_timestamp = int(time.time())
	local_time = str(int(round(time.time() * 1000)))
	xx = basename + local_time + ".dat"
	if os.path.isfile(xx):
		time.sleep(1)
		return generate_filename(basename)
	return xx

#..................................................Traslations MEAN
def Trasl_Mean(A):
	a = A.shape
	B = np.zeros((a[1],a[1]), dtype=np.float)
	for i in range(a[1]):
		B[i] = np.roll(A[i],-i)
	return np.mean(B, axis=0)







