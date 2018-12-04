#!/bin/bash

### nom du job (a changer)
#$ -N N_N1N-L_L1L-U_U1U-O_O1O

### parallel environment & nb cpu (NSLOTS)

#$ -pe mpi16_debian 16

#$ -q   "	h48-E5-2667v2deb128,h6-E5-2667v4deb128,h6-E5-2667v4deb128,E5-2667v2deb128nl,E5-2667v4deb256A,E5-2670deb128C,E5-2670deb128D,E5-2670deb128A,E5-2670deb128B,E5-2670deb128nl"

module load Python/3.6.1

### exporter les variables d'environnement sur tous les noeuds d'execution
#$ -V

###$ -e /dev/null
###$ -o /dev/null
 
WORKDIR="/home/pnaldesi/exact_di/ex_di/code"
cd ${WORKDIR}

echo $SHELL

mpirun -np 30 python3 bose.py N1N L1L U1U O1O

wait

rm -f uga-N_N1N-L_L1L-U_U1U-O_O1O.inp