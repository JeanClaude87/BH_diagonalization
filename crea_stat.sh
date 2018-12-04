rm -f script_wait.sh
rm -f u*

python3 setup.py build_ext --inplace
wait    

printf "%s\n" "#!/bin/bash" >> script_wait.sh

for NN in 2 3 4 5
    do
    for LL in 10 15 20 25
        do
        for UU in $(seq -w 0.0 0.2 1.0)
            do
            for OO in $(seq -w 0.0 0.02 0.5)
                do

                    sed -e "s/N1N/$NN/g" -e "s/L1L/$LL/g" -e "s/U1U/$UU/g" -e "s/O1O/$OO/g" < lancio.sh > temp.tmp

                    mv temp.tmp uga-N_$NN-L_$LL-U_$UU-O_$OO.inp        
                    printf  "%s\n" "uga-N_$NN-L_$LL-U_$UU-O_$OO.inp" >> script_wait.sh
                    
                    wait                

                done
            done
        done
    done


chmod +x script_wait.sh
#./script_wait.sh



