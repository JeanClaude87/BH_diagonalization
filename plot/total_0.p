set autoscale


	set terminal postscript fontfile add '/Users/naldesi/Type1/sfrm1200.pfb' \
							fontfile add '/Users/naldesi/Type1/sfti1200.pfb' \
							fontfile add '/Users/naldesi/Type1/cmmi10.pfb'
							
							
	set terminal postscript	eps color enhanced font "SFRM1200" 22	size 6,7
	set border lw 2 front

	set output	'total_0.eps'		
 	
########   Fitting function    ######################################################################################

 first(x) = ($0 > 0 ? base : base = x)
 first1(x) = ($0 > 0 ? base : base = abs(x))

	set multiplot layout 3,2

    set style line 1	lt 1 	lw 3 pt 2 	ps 1.4 lc rgb '#a6cee3' 
    set style line 2	lt 2	lw 3 pt 4 	ps 1.4 lc rgb '#1f78b4' 
    set style line 3	lt 3 	lw 3 pt 6 	ps 1.4 lc rgb '#b2df8a' 
    set style line 4	lt 4 	lw 3 pt 8 	ps 1.4 lc rgb '#33a02c' 
    set style line 5	lt 5 	lw 3 pt 10 	ps 1.4 lc rgb '#394bb8' 
    set style line 6	lt 6 	lw 3 pt 12 	ps 1.4 lc rgb '#c80236' 
    set style line 7	lt 7 	lw 3 pt 2 	ps 1.4 lc rgb '#fdbf6f'
    set style line 8	lt 8 	lw 3 pt 4 	ps 1.4 lc rgb '#ff7f00' 
	set style line 9	lt 9 	lw 3 pt 6 	ps 1.4 lc rgb '#cab2d6' 
	set style line 10	lt 10 	lw 3 pt 8 	ps 1.4 lc rgb '#6a3d9a' 
	set style line 11	lt 11 	lw 3 pt 10 	ps 1.4 lc rgb '#b15928' 
	set style line 16	lt 1 	lw 3 pt 12 	ps 1.4 lc rgb '#e31a1c' 

	set xlabel 't/J'	offset 0.0,0.0	norotate 	font "SFTI1200,28" 
	set xtic 2000	offset 0.0,0.0	 mirror		font "SFRM1200,22"
	set mxtic			

	set ylabel ''	offset 0.0,0.0 	rotate 		font "SFTI1200,28" 
	set ytic		offset 0.0,0.0	 mirror		font "SFRM1200,22"
	set mytic	

	set xrange[*:10000]
	set yrange[0:2.2]


	set title 'N=2 L=10 U/J=-5 V_0=0.007'
	plot 	'../dati_0/L_10/N_2/U_-5.0/bb_0.007/fidelity_cat_s.dat'	w l ls 4 t 'fid',\
			'../dati_0/L_10/N_2/U_-5.0/bb_0.007/corrente.dat'			w l ls 6 t 'cur',\
			1.0 w l dt 1 lc 'black' notitle,\
			0.5 w l dt 3 lc 'black' notitle


	set title 'N=3 L=10 U/J=-3 V_0=0.003'
	plot 	'../dati_0/L_10/N_3/U_-3.0/bb_0.003/fidelity_cat_s.dat'	w l ls 4 t 'fid',\
			'../dati_0/L_10/N_3/U_-3.0/bb_0.003/corrente.dat'			w l ls 6 t 'cur',\
			1.0 w l dt 1 lc 'black' notitle,\
			0.5 w l dt 3 lc 'black' notitle


	set title 'N=4 L=10 U/J=-2 V_0=0.001'
	plot 	'../dati_0/L_10/N_4/U_-2.0/bb_0.0007/fidelity_cat_s.dat'	w l ls 4 t 'fid',\
			'../dati_0/L_10/N_4/U_-2.0/bb_0.0007/corrente.dat'		w l ls 6 t 'cur',\
			1.0 w l dt 1 lc 'black' notitle,\
			0.5 w l dt 3 lc 'black' notitle


	set title 'N=5 L=10 U/J=-1.5 V_0=0.001'
	plot 	'../dati_0/L_10/N_5/U_-1.5/bb_0.001/fidelity_cat_s.dat'	w l ls 4 t 'fid',\
			'../dati_0/L_10/N_5/U_-1.5/bb_0.001/corrente.dat'			w l ls 6 t 'cur',\
			1.0 w l dt 1 lc 'black' notitle,\
			0.5 w l dt 3 lc 'black' notitle


	set title 'N=6 L=10 U/J=-1.5 V_0=0.0007'
	plot 	'../dati_0/L_10/N_6/U_-1.0/bb_0.0007/fidelity_cat_s.dat'	w l ls 4 t 'fid',\
			'../dati_0/L_10/N_6/U_-1.0/bb_0.0007/corrente.dat'		w l ls 6 t 'cur',\
			1.0 w l dt 1 lc 'black' notitle,\
			0.5 w l dt 3 lc 'black' notitle

	unset multiplot
	set output	

	
	!pstopdf total_0.eps
	!rm total_0.eps


