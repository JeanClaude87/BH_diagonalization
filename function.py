import numpy as np
import itertools
import os
from scipy.sparse import csc_matrix
import time


#..................................hilbert space dimension
def fact_creation(**args):

	ll = args.get("ll")
	nn = args.get("nn")	

	tab = [1 for x in range(nn+ll+2)]

	for i in range(nn+ll-1):
		tab[i+1] = int(tab[i]*(i+1))

	return tab


def hilb_dim(x,y,tab_fact):
	
	#x = args.get("nn")
	#y = args.get("ll")	
	#tab_fact = args.get("tab_fact")

	kk  = tab_fact[x+y-1]/(tab_fact[x]*tab_fact[y-1])

	uga = int(kk)
	
	return uga

def hilb_dim_tab(**args):

	ll = args.get("ll")
	nn = args.get("nn")	
	tab_fact = args.get("tab_fact")

	tab = [[hilb_dim(n, l, tab_fact) for l in range(ll+1)] for n in range(nn+1)]

	return tab

#..................................base preparation

def Base_prep(**args):

	ll = args.get("ll")
	nn = args.get("nn")

	base_bin = []
	base_num = []
	base_ind = []

#	Max = int(sum([2**(ll+nn-2-x) for x in range(nn)]))
#	TO_con_tab = [None] * Max

	for bits in itertools.combinations(range(nn+ll-1), nn):
		s = ['0'] * (nn+ll-1)	
		for bit in bits:
			s[bit] = '1'	
		bi = ''.join(s)	
		base_bin.append(bi)

		bose = TO_bose_conf(s,ll)

		base_num.append(bose)

		in_bi=int(bi,2)
		base_ind.append(in_bi)

		#TO_con_tab[in_bi] = bi

	base_bose = np.asarray(base_num, dtype=np.int8)

	return base_bin, base_bose, base_ind


#..................................index search
def get_index(state,**args):

	ll           = args.get("ll")
	nn           = args.get("nn")
	DIM_H        = args.get("DIM_H")
	tab_fact     = args.get("tab_fact")
	hilb_dim_tab = args.get("hilb_dim_tab")

	size   = int(nn+ll-1)
	r_par  = int(nn)    #remaining_particles
	r_sit  = int(ll)	#remaining_sites
	result = DIM_H

	for jj in range(size):
		action_i = int(state[jj]);

		if r_par==0:
			break
		if action_i == 0:
			#print(jj,r_par)		
		
			result -= hilb_dim_tab[r_par-1][r_sit]
			#result -= binomial_table[r_par-1][r_sit]
			r_sit  -= 1;

		else:
			r_par -= 1;
			#print('else',jj,r_par)

	return DIM_H-result


#..................................from BOSE configuration to bin number
def FROM_bose_TO_bin(state,**args):

	i=[]

	for x in state:
		if x == 0:
			i.append('0')
		else:
			for x in range(x):
				i.append('1')
			i.append('0')

	i   = i[:-1]
	uga = ''.join(i)

	return uga


#..................................from configuration to bin number
def TO_bin(xx):
	return int(xx,2)

#..................................from bin number to configuration
def TO_con(x,L):
	x1=int(x)
	L1=int(L)
	bi=np.binary_repr(x1, width=L1)

#	print(x1,bi)

	return  bi

#def TO_con_1(x,L):
	
#	ind = CONF_tab_sp[x]

	#ind0 = CONF_tab_sp[x].toarray()
	#ind  = ind0[0,0]

#	return BASE_bin[ind]

#..................................hop. preparation
#............. ......BC=0 periodic, BC=1 open
def Hop_prep(**args):

	X = args.get("ll")
	Y = args.get("nn")

	Hop_dim=X+Y-2
	
	return [TO_con(2**i+2**((i+1)%(X+Y-1)),X+Y-1) for i in range(Hop_dim)]

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


#..................................................generate filename
def generate_filename(basename):
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


#..................................................Print_MATRIX
def print_matrix(H):

	#print('matrix to print')

	if isinstance(H, csc_matrix):
		print_h = csc_matrix.todense(H)
		print(print_h)
	else:
		print(H)

	return 0














