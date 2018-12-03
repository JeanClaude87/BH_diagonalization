#!/bin/bash

for NN in 2 3
	do
	for LL in 5 10
		do
		for UU in $(seq -w 0.0 0.1 0.3)
			do
			for OO in $(seq -w 0.0 0.1 0.3)
				do

					mpirun -np 1 python3 bose.py $NN $LL $UU $OO
					wait

			done
		done
	done
done

