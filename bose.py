import numpy as np
from math import factorial
import math
import itertools

import function as ff


ll=10
nn=5

tab_fact = ff.fact_creation(nn+ll)
DIM_H    = ff.hilb_dim(tab_fact,nn,ll)

BC=0

BASE_num, BASE_bin = ff.Base_prep(ll,nn)

HAM = ff.bose_sparse_Hamiltonian(ll,nn,BC,BASE_bin,tab_fact)

