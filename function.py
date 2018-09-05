import numpy as np
from math import factorial
import math
import itertools


#..................................hilbert space dimension
def fact_creation(a):
	tab = [factorial(x) for x in range(a)]
	return tab

def fast_fact(fact_tab,x):
	return fact_tab[x]


def hilb_dim(fact_tab,n,l):
	kk = fast_fact(fact_tab,n+l-1)/(fast_fact(fact_tab,n)*fast_fact(fact_tab,l-1)) 
	uga= int(kk)
	return uga

#..................................base preparation
# n number of particles
# l number of sites
def Base_prep(l,n):
	base_bin = []
	base_num = []
	for bits in itertools.combinations(range(n+l-1), n):
		s = ['0'] * (n+l-1)	
		for bit in bits:
			s[bit] = '1'
		base_num.append(s)
		base_bin.append(''.join(s))
	return base_num, base_bin

#..................................index search
# n number of particles
# l number of sites
def get_index(state, l, n, fact_tab):

	size = int(n+l-1)
	r_par = int(n)  #remaining_particles
	r_sit = int(l)	#remaining_sites
	result = hilb_dim(fact_tab,n,l)

	#print(size, len(state), state)

	for jj in range(size):
		action_i = int(state[jj]);
		if r_par==0:
			break
		if action_i == 0:
			#print(jj,r_par)		
		
			result -= hilb_dim(fact_tab, r_par-1, r_sit);
			r_sit  -= 1;

		else:
			r_par -= 1;
			#print('else',jj,r_par)

	return hilb_dim(fact_tab,n,l)-result

#..................................from configuration to bin number
def TO_bin(xx):
	return int(xx,2)

#..................................from bin number to configuration
def TO_con(x,L):
	x1=int(x)
	L1=int(L)
	return np.binary_repr(x1, width=L1)

#..................................hop. preparation
#...................BC=0 periodic, BC=1 open
def Hop_prep(L,N):
	Hop_dim=L+N-2
	return [TO_con(2**i+2**((i+1)%(L+N-1)),L+N-1) for i in range(Hop_dim)]

#..................................counting number of one
POPCOUNT_TABLE16 = [0] * 2**16
for index in range(len(POPCOUNT_TABLE16)):
	POPCOUNT_TABLE16[index] = (index & 1) + POPCOUNT_TABLE16[index >> 1]

def one_count(v):
	return (POPCOUNT_TABLE16[ v        & 0xffff] +
			POPCOUNT_TABLE16[(v >> 16) & 0xffff])

def TO_bose_conf(x,ll):
	p=0
	conf = np.zeros(ll, dtype=np.int)
	for jj in range(len(x)):
		if x[jj]=='1':
			conf[p]+=1
		else: 
			p+=1	
	return 	conf


def kin_expval(x,y,ll,BASE_bin):

	state_x   = BASE_bin[x]
	bosecon_x = TO_bose_conf(state_x,ll)

	state_y   = BASE_bin[y]
	bosecon_y = TO_bose_conf(state_y,ll)

	jump_c  =	np.argmax(bosecon_x-bosecon_y)
	jump_cd  =	np.argmin(bosecon_x-bosecon_y)

	result	  = np.sqrt((bosecon_x[jump_cd]+1)*(bosecon_x[jump_c])) 
	
	return result


def bose_Hamiltonian (ll,nn,BC,t,U,BASE_bin,tab_fact):

	DIM_H=hilb_dim(tab_fact,nn,ll)
	ham_ind1 = []
	ham_ind2 = []
	ham_val  = []

	for i in range(DIM_H):
		state = BASE_bin[i]

##----- INTERACTIONS

		bosecon = TO_bose_conf(state,ll)

		int_val = 0
		for x in range(ll):
			nx = bosecon[x]
			int_val += U*nx*(nx-1.)/(2.)

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

				j = get_index(hop_state , ll, nn, tab_fact)	
				kin_val = t*kin_expval(i,j,ll,BASE_bin)

		#---- KINETIC = we store i,j,kin_val !!!!
				ham_ind1.append( i )
				ham_ind2.append( j )
				ham_val.append( kin_val )
	
		# here we put the PERIODIC
		#...................BC=0 periodic, BC=1 open

		if BC == 0:
			if state[0] == '1':

				PBC_newstate1 = state[1:]+state[0]				
				j = get_index(PBC_newstate1, ll, nn, tab_fact)	

				kin_val = t*kin_expval(i,j,ll,BASE_bin)
				
		#---- KINETIC = we store i,j,kin_val !!!!
				ham_ind1.append( i )
				ham_ind2.append( j )
				ham_val.append( kin_val )

			if state[-1] == '1':

				PBC_newstate2 = state[-1]+state[:-1]
				j = get_index(PBC_newstate2, ll, nn, tab_fact)	

				kin_val = t*kin_expval(i,j,ll,BASE_bin)
				
		#---- KINETIC = we store i,j,kin_val !!!!
				ham_ind1.append( i )
				ham_ind2.append( j )
				ham_val.append( kin_val )
		
	return ham_ind1,ham_ind2,ham_val
















