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
def get_index(state, l, n):

	size = int(n+l-1)
	r_par = int(n)  #remaining_particles
	r_sit = int(l)	#remaining_sites
	result = hilb_dim(n,l)

	for jj in range(size):
		action_i = int(state[jj]);
		if r_par==0:
			break
		if action_i == 0:
			#print(jj,r_par)		
		
			result -= hilb_dim(r_par-1,r_sit);
			r_sit  -= 1;

		else:
			r_par -= 1;
			#print('else',jj,r_par)

	return result 

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
def Hop_prep(L,N,BC):
	if BC == 1:
		Hop_dim=L+N-1
	else:
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

def bose_sparse_Hamiltonian (ll,nn,BC,BASE_bin,tab_fact):

	for hh in range(hilb_dim(tab_fact,nn,ll)):
		state = BASE_bin[hh]
		#
		# -> we fill the INTERACTION
		#
		for hop in Hop_prep(ll,nn,BC):
			uga = TO_bin(state)^TO_bin(hop)
			if one_count(uga) == nn:
				ciao=0
				#print(TO_bose_conf(state,ll),
				#	hop,
				#	TO_bose_conf(TO_con(uga,ll+nn-1),ll))

				#
				# -> we fill the KINETIC
				#
	return 0
