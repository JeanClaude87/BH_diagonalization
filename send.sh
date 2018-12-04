#!/bin/bash

for NN in 4 5
	do
	for LL in 10 15 20 25
		do
		for UU in $(seq -w 0.0 0.2 1.0)
			do
			for OO in $(seq -w 0.0 0.02 0.5)
				do

				mpirun -np 30 python3 bose.py $NN $LL $UU $OO
				echo $NN $LL $UU $OO 


			done
		done
	done
done
