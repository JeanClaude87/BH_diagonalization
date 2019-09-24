set autoscale


	set terminal postscript fontfile add '/Users/naldesi/Type1/sfrm1200.pfb' \
							fontfile add '/Users/naldesi/Type1/sfti1200.pfb' \
							fontfile add '/Users/naldesi/Type1/cmmi10.pfb'
							
							
	set terminal postscript	eps color enhanced font "SFRM1200" 22	size 18,20
	set border lw 2 front

	set output	'total.eps'		
 	
#######   Fitting function    ######################################################################################

 first(x) = ($0 > 0 ? base : base = x)
 first1(x) = ($0 > 0 ? base : base = abs(x))

	set multiplot layout 4,3

	set key font "SFRM1200,22" vertical maxrows 3 samplen 1.6 spacing 1.2

    set style line 1	lt 1 	lw 3 pt 2 	ps 1.0 lc rgb '#a6cee3' 
    set style line 2	lt 2	lw 3 pt 4 	ps 1.0 lc rgb '#1f78b4' 
    set style line 3	lt 3 	lw 3 pt 6 	ps 1.0 lc rgb '#b2df8a' 
    set style line 4	lt 4 	lw 3 pt 8 	ps 1.0 lc rgb '#33a02c' 
    set style line 5	lt 5 	lw 3 pt 10 	ps 1.0 lc rgb '#394bb8' 
    set style line 6	lt 6 	lw 3 pt 12 	ps 1.0 lc rgb '#c80236' 
    set style line 7	lt 7 	lw 3 pt 2 	ps 1.0 lc rgb '#fdbf6f'
    set style line 8	lt 8 	lw 3 pt 4 	ps 1.0 lc rgb '#ff7f00' 
	set style line 9	lt 9 	lw 3 pt 6 	ps 1.0 lc rgb '#cab2d6' 
	set style line 10	lt 10 	lw 3 pt 8 	ps 1.0 lc rgb '#6a3d9a' 
	set style line 11	lt 11 	lw 3 pt 10 	ps 1.0 lc rgb '#b15928' 
	set style line 16	lt 1 	lw 3 pt 12 	ps 1.0 lc rgb '#e31a1c' 

	set xlabel 't/J'	offset 0.0,0.0	norotate 	font "SFTI1200,28" 
	set xtic 1000	offset 0.0,0.0	 mirror		font "SFRM1200,22"
	set mxtic			

	set ylabel 'fidelity'	offset 0.0,0.0 	rotate 		font "SFTI1200,28" 
	set ytic		offset 0.0,0.0	 mirror		font "SFRM1200,22"
	set mytic	

	set xrange[*:5000]
	set yrange[*:1.5]

	set title 'N=2 L=10 U/J=-5 SYMMETRIC CAT'
	plot 	'../dati/L_10/N_2/U_-5.0/bb_0.05/fidelity_cat_s.dat' 	w l ls 2 t 'V_0 = 0.05' ,\
			'../dati/L_10/N_2/U_-5.0/bb_0.03/fidelity_cat_s.dat' 	w l ls 3 t 'V_0 = 0.03' ,\
			'../dati/L_10/N_2/U_-5.0/bb_0.01/fidelity_cat_s.dat' 	w l ls 4 t 'V_0 = 0.01' ,\
			'../dati/L_10/N_2/U_-5.0/bb_0.007/fidelity_cat_s.dat' w l ls 5 t 'V_0 = 0.007' ,\
			'../dati/L_10/N_2/U_-5.0/bb_0.005/fidelity_cat_s.dat' w l ls 6 t 'V_0 = 0.005' ,\
			'../dati/L_10/N_2/U_-5.0/bb_0.003/fidelity_cat_s.dat' w l ls 7 t 'V_0 = 0.003' ,\
			'../dati/L_10/N_2/U_-5.0/bb_0.001/fidelity_cat_s.dat' w l ls 8 t 'V_0 = 0.001' ,\
			'../dati/L_10/N_2/U_-5.0/bb_0.0007/fidelity_cat_s.dat' w l ls 9 t 'V_0 = 0.0007' ,\
			'../dati/L_10/N_2/U_-5.0/bb_0.0005/fidelity_cat_s.dat' w l ls 1 t 'V_0 = 0.0005' ,\
			'../dati/L_10/N_2/U_-5.0/bb_0.0003/fidelity_cat_s.dat' w l ls 2 t 'V_0 = 0.0003' ,\
			'../dati/L_10/N_2/U_-5.0/bb_0.0001/fidelity_cat_s.dat' w l ls 3 t 'V_0 = 0.0001' ,\
			1 w l ls 9 notitle

	set title 'N=2 L=10 U/J=-5 ANTISYMMETRIC CAT'
	plot 	'../dati/L_10/N_2/U_-5.0/bb_0.05/fidelity_cat_a.dat' 	w l ls 2 t 'V_0 = 0.05' ,\
			'../dati/L_10/N_2/U_-5.0/bb_0.03/fidelity_cat_a.dat' 	w l ls 3 t 'V_0 = 0.03' ,\
			'../dati/L_10/N_2/U_-5.0/bb_0.01/fidelity_cat_a.dat' 	w l ls 4 t 'V_0 = 0.01' ,\
			'../dati/L_10/N_2/U_-5.0/bb_0.007/fidelity_cat_a.dat' w l ls 5 t 'V_0 = 0.007' ,\
			'../dati/L_10/N_2/U_-5.0/bb_0.005/fidelity_cat_a.dat' w l ls 6 t 'V_0 = 0.005' ,\
			'../dati/L_10/N_2/U_-5.0/bb_0.003/fidelity_cat_a.dat' w l ls 7 t 'V_0 = 0.003' ,\
			'../dati/L_10/N_2/U_-5.0/bb_0.001/fidelity_cat_a.dat' w l ls 8 t 'V_0 = 0.001' ,\
			'../dati/L_10/N_2/U_-5.0/bb_0.0007/fidelity_cat_a.dat' w l ls 9 t 'V_0 = 0.0007' ,\
			'../dati/L_10/N_2/U_-5.0/bb_0.0005/fidelity_cat_a.dat' w l ls 1 t 'V_0 = 0.0005' ,\
			'../dati/L_10/N_2/U_-5.0/bb_0.0003/fidelity_cat_a.dat' w l ls 2 t 'V_0 = 0.0003' ,\
			'../dati/L_10/N_2/U_-5.0/bb_0.0001/fidelity_cat_a.dat' w l ls 3 t 'V_0 = 0.0001' ,\
			1 w l ls 9 notitle 

	set ylabel 'current'	offset 0.0,0.0 	rotate 		font "SFTI1200,28" 
	set ytic				offset 0.0,0.0	 mirror		font "SFRM1200,22"
	set mytic	

	set yrange[*:3]

	set title 'N=2 L=10 U/J=-5 CURRENT'
	plot 	'../dati/L_10/N_2/U_-5.0/bb_0.05/corrente.dat' 	w l ls 2 t 'V_0 = 0.05' ,\
			'../dati/L_10/N_2/U_-5.0/bb_0.03/corrente.dat' 	w l ls 3 t 'V_0 = 0.03' ,\
			'../dati/L_10/N_2/U_-5.0/bb_0.01/corrente.dat' 	w l ls 4 t 'V_0 = 0.01' ,\
			'../dati/L_10/N_2/U_-5.0/bb_0.007/corrente.dat' w l ls 5 t 'V_0 = 0.007' ,\
			'../dati/L_10/N_2/U_-5.0/bb_0.005/corrente.dat' w l ls 6 t 'V_0 = 0.005' ,\
			'../dati/L_10/N_2/U_-5.0/bb_0.003/corrente.dat' w l ls 7 t 'V_0 = 0.003' ,\
			'../dati/L_10/N_2/U_-5.0/bb_0.001/corrente.dat' w l ls 8 t 'V_0 = 0.001' ,\
			'../dati/L_10/N_2/U_-5.0/bb_0.0007/corrente.dat' w l ls 9 t 'V_0 = 0.0007' ,\
			'../dati/L_10/N_2/U_-5.0/bb_0.0005/corrente.dat' w l ls 1 t 'V_0 = 0.0005' ,\
			'../dati/L_10/N_2/U_-5.0/bb_0.0003/corrente.dat' w l ls 2 t 'V_0 = 0.0003' ,\
			'../dati/L_10/N_2/U_-5.0/bb_0.0001/corrente.dat' w l ls 3 t 'V_0 = 0.0001' ,\
			1 w l ls 9 notitle 


	set ylabel 'fidelity'	offset 0.0,0.0 	rotate 		font "SFTI1200,28" 
	set ytic		offset 0.0,0.0	 mirror		font "SFRM1200,22"
	set mytic	

	set yrange[*:1.5]

	set title 'N=3 L=10 U/J=-3 SYMMETRIC CAT'
	plot 	'../dati/L_10/N_3/U_-3.0/bb_0.05/fidelity_cat_s.dat' 	w l ls 2 t 'V_0 = 0.05' ,\
			'../dati/L_10/N_3/U_-3.0/bb_0.03/fidelity_cat_s.dat' 	w l ls 3 t 'V_0 = 0.03' ,\
			'../dati/L_10/N_3/U_-3.0/bb_0.01/fidelity_cat_s.dat' 	w l ls 4 t 'V_0 = 0.01' ,\
			'../dati/L_10/N_3/U_-3.0/bb_0.007/fidelity_cat_s.dat' w l ls 5 t 'V_0 = 0.007' ,\
			'../dati/L_10/N_3/U_-3.0/bb_0.005/fidelity_cat_s.dat' w l ls 6 t 'V_0 = 0.005' ,\
			'../dati/L_10/N_3/U_-3.0/bb_0.003/fidelity_cat_s.dat' w l ls 7 t 'V_0 = 0.003' ,\
			'../dati/L_10/N_3/U_-3.0/bb_0.001/fidelity_cat_s.dat' w l ls 8 t 'V_0 = 0.001' ,\
			'../dati/L_10/N_3/U_-3.0/bb_0.0007/fidelity_cat_s.dat' w l ls 9 t 'V_0 = 0.0007' ,\
			'../dati/L_10/N_3/U_-3.0/bb_0.0005/fidelity_cat_s.dat' w l ls 1 t 'V_0 = 0.0005' ,\
			'../dati/L_10/N_3/U_-3.0/bb_0.0003/fidelity_cat_s.dat' w l ls 2 t 'V_0 = 0.0003' ,\
			'../dati/L_10/N_3/U_-3.0/bb_0.0001/fidelity_cat_s.dat' w l ls 3 t 'V_0 = 0.0001' ,\
			1 w l ls 9 notitle 

	set title 'N=3 L=10 U/J=-3 ANTISYMMETRIC CAT'
	plot 	'../dati/L_10/N_3/U_-3.0/bb_0.05/fidelity_cat_a.dat' 	w l ls 2 t 'V_0 = 0.05' ,\
			'../dati/L_10/N_3/U_-3.0/bb_0.03/fidelity_cat_a.dat' 	w l ls 3 t 'V_0 = 0.03' ,\
			'../dati/L_10/N_3/U_-3.0/bb_0.01/fidelity_cat_a.dat' 	w l ls 4 t 'V_0 = 0.01' ,\
			'../dati/L_10/N_3/U_-3.0/bb_0.007/fidelity_cat_a.dat' w l ls 5 t 'V_0 = 0.007' ,\
			'../dati/L_10/N_3/U_-3.0/bb_0.005/fidelity_cat_a.dat' w l ls 6 t 'V_0 = 0.005' ,\
			'../dati/L_10/N_3/U_-3.0/bb_0.003/fidelity_cat_a.dat' w l ls 7 t 'V_0 = 0.003' ,\
			'../dati/L_10/N_3/U_-3.0/bb_0.001/fidelity_cat_a.dat' w l ls 8 t 'V_0 = 0.001' ,\
			'../dati/L_10/N_3/U_-3.0/bb_0.0007/fidelity_cat_a.dat' w l ls 9 t 'V_0 = 0.0007' ,\
			'../dati/L_10/N_3/U_-3.0/bb_0.0005/fidelity_cat_a.dat' w l ls 1 t 'V_0 = 0.0005' ,\
			'../dati/L_10/N_3/U_-3.0/bb_0.0003/fidelity_cat_a.dat' w l ls 2 t 'V_0 = 0.0003' ,\
			'../dati/L_10/N_3/U_-3.0/bb_0.0001/fidelity_cat_a.dat' w l ls 3 t 'V_0 = 0.0001' ,\
			1 w l ls 9 notitle 

	set ylabel 'current'	offset 0.0,0.0 	rotate 		font "SFTI1200,28" 
	set ytic				offset 0.0,0.0	 mirror		font "SFRM1200,22"
	set mytic	

	set yrange[*:3]

	set title 'N=3 L=10 U/J=-3 CURRENT'
	plot 	'../dati/L_10/N_3/U_-3.0/bb_0.05/corrente.dat' 	w l ls 2 t 'V_0 = 0.05' ,\
			'../dati/L_10/N_3/U_-3.0/bb_0.03/corrente.dat' 	w l ls 3 t 'V_0 = 0.03' ,\
			'../dati/L_10/N_3/U_-3.0/bb_0.01/corrente.dat' 	w l ls 4 t 'V_0 = 0.01' ,\
			'../dati/L_10/N_3/U_-3.0/bb_0.007/corrente.dat' w l ls 5 t 'V_0 = 0.007' ,\
			'../dati/L_10/N_3/U_-3.0/bb_0.005/corrente.dat' w l ls 6 t 'V_0 = 0.005' ,\
			'../dati/L_10/N_3/U_-3.0/bb_0.003/corrente.dat' w l ls 7 t 'V_0 = 0.003' ,\
			'../dati/L_10/N_3/U_-3.0/bb_0.001/corrente.dat' w l ls 8 t 'V_0 = 0.001' ,\
			'../dati/L_10/N_3/U_-3.0/bb_0.0007/corrente.dat' w l ls 9 t 'V_0 = 0.0007' ,\
			'../dati/L_10/N_3/U_-3.0/bb_0.0005/corrente.dat' w l ls 1 t 'V_0 = 0.0005' ,\
			'../dati/L_10/N_3/U_-3.0/bb_0.0003/corrente.dat' w l ls 2 t 'V_0 = 0.0003' ,\
			'../dati/L_10/N_3/U_-3.0/bb_0.0001/corrente.dat' w l ls 3 t 'V_0 = 0.0001' ,\
			1 w l ls 9 notitle 


	set ylabel 'fidelity'	offset 0.0,0.0 	rotate 		font "SFTI1200,28" 
	set ytic		offset 0.0,0.0	 mirror		font "SFRM1200,22"
	set mytic	

	set yrange[*:1.5]

	set title 'N=4 L=10 U/J=-2 SYMMETRIC CAT'
	plot 	'../dati/L_10/N_4/U_-2.0/bb_0.05/fidelity_cat_s.dat' 	w l ls 2 t 'V_0 = 0.05' ,\
			'../dati/L_10/N_4/U_-2.0/bb_0.03/fidelity_cat_s.dat' 	w l ls 3 t 'V_0 = 0.03' ,\
			'../dati/L_10/N_4/U_-2.0/bb_0.01/fidelity_cat_s.dat' 	w l ls 4 t 'V_0 = 0.01' ,\
			'../dati/L_10/N_4/U_-2.0/bb_0.007/fidelity_cat_s.dat' w l ls 5 t 'V_0 = 0.007' ,\
			'../dati/L_10/N_4/U_-2.0/bb_0.005/fidelity_cat_s.dat' w l ls 6 t 'V_0 = 0.005' ,\
			'../dati/L_10/N_4/U_-2.0/bb_0.003/fidelity_cat_s.dat' w l ls 7 t 'V_0 = 0.003' ,\
			'../dati/L_10/N_4/U_-2.0/bb_0.001/fidelity_cat_s.dat' w l ls 8 t 'V_0 = 0.001' ,\
			'../dati/L_10/N_4/U_-2.0/bb_0.0007/fidelity_cat_s.dat' w l ls 9 t 'V_0 = 0.0007' ,\
			'../dati/L_10/N_4/U_-2.0/bb_0.0005/fidelity_cat_s.dat' w l ls 1 t 'V_0 = 0.0005' ,\
			'../dati/L_10/N_4/U_-2.0/bb_0.0003/fidelity_cat_s.dat' w l ls 2 t 'V_0 = 0.0003' ,\
			'../dati/L_10/N_4/U_-2.0/bb_0.0001/fidelity_cat_s.dat' w l ls 3 t 'V_0 = 0.0001' ,\
			1 w l ls 9 notitle

	set title 'N=4 L=10 U/J=-2 ANTISYMMETRIC CAT'
	plot 	'../dati/L_10/N_4/U_-2.0/bb_0.05/fidelity_cat_a.dat' 	w l ls 2 t 'V_0 = 0.05' ,\
			'../dati/L_10/N_4/U_-2.0/bb_0.03/fidelity_cat_a.dat' 	w l ls 3 t 'V_0 = 0.03' ,\
			'../dati/L_10/N_4/U_-2.0/bb_0.01/fidelity_cat_a.dat' 	w l ls 4 t 'V_0 = 0.01' ,\
			'../dati/L_10/N_4/U_-2.0/bb_0.007/fidelity_cat_a.dat' w l ls 5 t 'V_0 = 0.007' ,\
			'../dati/L_10/N_4/U_-2.0/bb_0.005/fidelity_cat_a.dat' w l ls 6 t 'V_0 = 0.005' ,\
			'../dati/L_10/N_4/U_-2.0/bb_0.003/fidelity_cat_a.dat' w l ls 7 t 'V_0 = 0.003' ,\
			'../dati/L_10/N_4/U_-2.0/bb_0.001/fidelity_cat_a.dat' w l ls 8 t 'V_0 = 0.001' ,\
			'../dati/L_10/N_4/U_-2.0/bb_0.0007/fidelity_cat_a.dat' w l ls 9 t 'V_0 = 0.0007' ,\
			'../dati/L_10/N_4/U_-2.0/bb_0.0005/fidelity_cat_a.dat' w l ls 1 t 'V_0 = 0.0005' ,\
			'../dati/L_10/N_4/U_-2.0/bb_0.0003/fidelity_cat_a.dat' w l ls 2 t 'V_0 = 0.0003' ,\
			'../dati/L_10/N_4/U_-2.0/bb_0.0001/fidelity_cat_a.dat' w l ls 3 t 'V_0 = 0.0001' ,\
			1 w l ls 9 notitle

	set ylabel 'current'	offset 0.0,0.0 	rotate 		font "SFTI1200,28" 
	set ytic				offset 0.0,0.0	 mirror		font "SFRM1200,22"
	set mytic	

	set yrange[*:3]

	set title 'N=4 L=10 U/J=-2 CURRENT'
	plot 	'../dati/L_10/N_4/U_-2.0/bb_0.05/corrente.dat' 	w l ls 2 t 'V_0 = 0.05' ,\
			'../dati/L_10/N_4/U_-2.0/bb_0.03/corrente.dat' 	w l ls 3 t 'V_0 = 0.03' ,\
			'../dati/L_10/N_4/U_-2.0/bb_0.01/corrente.dat' 	w l ls 4 t 'V_0 = 0.01' ,\
			'../dati/L_10/N_4/U_-2.0/bb_0.007/corrente.dat' w l ls 5 t 'V_0 = 0.007' ,\
			'../dati/L_10/N_4/U_-2.0/bb_0.005/corrente.dat' w l ls 6 t 'V_0 = 0.005' ,\
			'../dati/L_10/N_4/U_-2.0/bb_0.003/corrente.dat' w l ls 7 t 'V_0 = 0.003' ,\
			'../dati/L_10/N_4/U_-2.0/bb_0.001/corrente.dat' w l ls 8 t 'V_0 = 0.001' ,\
			'../dati/L_10/N_4/U_-2.0/bb_0.0007/corrente.dat' w l ls 9 t 'V_0 = 0.0007' ,\
			'../dati/L_10/N_4/U_-2.0/bb_0.0005/corrente.dat' w l ls 1 t 'V_0 = 0.0005' ,\
			'../dati/L_10/N_4/U_-2.0/bb_0.0003/corrente.dat' w l ls 2 t 'V_0 = 0.0003' ,\
			'../dati/L_10/N_4/U_-2.0/bb_0.0001/corrente.dat' w l ls 3 t 'V_0 = 0.0001' ,\
			1 w l ls 9 notitle


	set ylabel 'fidelity'	offset 0.0,0.0 	rotate 		font "SFTI1200,28" 
	set ytic		offset 0.0,0.0	 mirror		font "SFRM1200,22"
	set mytic	

	set xrange[*:7000]
	set yrange[*:1.5]

	set title 'N=5 L=10 U/J=-1.5 SYMMETRIC CAT'
	plot 	'../dati/L_10/N_5/U_-1.5/bb_0.05/fidelity_cat_s.dat' 	w l ls 2 t 'V_0 = 0.05' ,\
			'../dati/L_10/N_5/U_-1.5/bb_0.03/fidelity_cat_s.dat' 	w l ls 3 t 'V_0 = 0.03' ,\
			'../dati/L_10/N_5/U_-1.5/bb_0.01/fidelity_cat_s.dat' 	w l ls 4 t 'V_0 = 0.01' ,\
			'../dati/L_10/N_5/U_-1.5/bb_0.007/fidelity_cat_s.dat' w l ls 5 t 'V_0 = 0.007' ,\
			'../dati/L_10/N_5/U_-1.5/bb_0.005/fidelity_cat_s.dat' w l ls 6 t 'V_0 = 0.005' ,\
			'../dati/L_10/N_5/U_-1.5/bb_0.003/fidelity_cat_s.dat' w l ls 7 t 'V_0 = 0.003' ,\
			'../dati/L_10/N_5/U_-1.5/bb_0.001/fidelity_cat_s.dat' w l ls 8 t 'V_0 = 0.001' ,\
			'../dati/L_10/N_5/U_-1.5/bb_0.0007/fidelity_cat_s.dat' w l ls 9 t 'V_0 = 0.0007' ,\
			'../dati/L_10/N_5/U_-1.5/bb_0.0005/fidelity_cat_s.dat' w l ls 1 t 'V_0 = 0.0005' ,\
			'../dati/L_10/N_5/U_-1.5/bb_0.0003/fidelity_cat_s.dat' w l ls 2 t 'V_0 = 0.0003' ,\
			'../dati/L_10/N_5/U_-1.5/bb_0.0001/fidelity_cat_s.dat' w l ls 3 t 'V_0 = 0.0001' ,\
			1 w l ls 9 notitle

	set title 'N=5 L=10 U/J=-1.5 ANTISYMMETRIC CAT'
	plot 	'../dati/L_10/N_5/U_-1.5/bb_0.05/fidelity_cat_a.dat' 	w l ls 2 t 'V_0 = 0.05' ,\
			'../dati/L_10/N_5/U_-1.5/bb_0.03/fidelity_cat_a.dat' 	w l ls 3 t 'V_0 = 0.03' ,\
			'../dati/L_10/N_5/U_-1.5/bb_0.01/fidelity_cat_a.dat' 	w l ls 4 t 'V_0 = 0.01' ,\
			'../dati/L_10/N_5/U_-1.5/bb_0.007/fidelity_cat_a.dat' w l ls 5 t 'V_0 = 0.007' ,\
			'../dati/L_10/N_5/U_-1.5/bb_0.005/fidelity_cat_a.dat' w l ls 6 t 'V_0 = 0.005' ,\
			'../dati/L_10/N_5/U_-1.5/bb_0.003/fidelity_cat_a.dat' w l ls 7 t 'V_0 = 0.003' ,\
			'../dati/L_10/N_5/U_-1.5/bb_0.001/fidelity_cat_a.dat' w l ls 8 t 'V_0 = 0.001' ,\
			'../dati/L_10/N_5/U_-1.5/bb_0.0007/fidelity_cat_a.dat' w l ls 9 t 'V_0 = 0.0007' ,\
			'../dati/L_10/N_5/U_-1.5/bb_0.0005/fidelity_cat_a.dat' w l ls 1 t 'V_0 = 0.0005' ,\
			'../dati/L_10/N_5/U_-1.5/bb_0.0003/fidelity_cat_a.dat' w l ls 2 t 'V_0 = 0.0003' ,\
			'../dati/L_10/N_5/U_-1.5/bb_0.0001/fidelity_cat_a.dat' w l ls 3 t 'V_0 = 0.0001' ,\
			1 w l ls 9 notitle

	set ylabel 'current'	offset 0.0,0.0 	rotate 		font "SFTI1200,28" 
	set ytic				offset 0.0,0.0	 mirror		font "SFRM1200,22"
	set mytic	

	set yrange[*:3]

	set title 'N=5 L=10 U/J=-1.5 CURRENT'
	plot 	'../dati/L_10/N_5/U_-1.5/bb_0.05/corrente.dat' 	w l ls 2 t 'V_0 = 0.05' ,\
			'../dati/L_10/N_5/U_-1.5/bb_0.03/corrente.dat' 	w l ls 3 t 'V_0 = 0.03' ,\
			'../dati/L_10/N_5/U_-1.5/bb_0.01/corrente.dat' 	w l ls 4 t 'V_0 = 0.01' ,\
			'../dati/L_10/N_5/U_-1.5/bb_0.007/corrente.dat' w l ls 5 t 'V_0 = 0.007' ,\
			'../dati/L_10/N_5/U_-1.5/bb_0.005/corrente.dat' w l ls 6 t 'V_0 = 0.005' ,\
			'../dati/L_10/N_5/U_-1.5/bb_0.003/corrente.dat' w l ls 7 t 'V_0 = 0.003' ,\
			'../dati/L_10/N_5/U_-1.5/bb_0.001/corrente.dat' w l ls 8 t 'V_0 = 0.001' ,\
			'../dati/L_10/N_5/U_-1.5/bb_0.0007/corrente.dat' w l ls 9 t 'V_0 = 0.0007' ,\
			'../dati/L_10/N_5/U_-1.5/bb_0.0005/corrente.dat' w l ls 1 t 'V_0 = 0.0005' ,\
			'../dati/L_10/N_5/U_-1.5/bb_0.0003/corrente.dat' w l ls 2 t 'V_0 = 0.0003' ,\
			'../dati/L_10/N_5/U_-1.5/bb_0.0001/corrente.dat' w l ls 3 t 'V_0 = 0.0001' ,\
			1 w l ls 9 notitle



	set ylabel 'fidelity'	offset 0.0,0.0 	rotate 		font "SFTI1200,28" 
	set ytic		offset 0.0,0.0	 mirror		font "SFRM1200,22"
	set mytic	

	set xrange[*:7000]
	set yrange[*:1.5]

#	set title 'N=6 L=10 U/J=-1 SYMMETRIC CAT'
#	plot 	'../dati/L_10/N_6/U_-1.0/bb_0.05/fidelity_cat_s.dat' 	w l ls 2 t 'V_0 = 0.05' ,\
			'../dati/L_10/N_6/U_-1.0/bb_0.03/fidelity_cat_s.dat' 	w l ls 3 t 'V_0 = 0.03' ,\
			'../dati/L_10/N_6/U_-1.0/bb_0.01/fidelity_cat_s.dat' 	w l ls 4 t 'V_0 = 0.01' ,\
			'../dati/L_10/N_6/U_-1.0/bb_0.007/fidelity_cat_s.dat' w l ls 5 t 'V_0 = 0.007' ,\
			'../dati/L_10/N_6/U_-1.0/bb_0.005/fidelity_cat_s.dat' w l ls 6 t 'V_0 = 0.005' ,\
			'../dati/L_10/N_6/U_-1.0/bb_0.003/fidelity_cat_s.dat' w l ls 7 t 'V_0 = 0.003' ,\
			'../dati/L_10/N_6/U_-1.0/bb_0.001/fidelity_cat_s.dat' w l ls 8 t 'V_0 = 0.001' ,\
			'../dati/L_10/N_6/U_-1.0/bb_0.0007/fidelity_cat_s.dat' w l ls 9 t 'V_0 = 0.0007' ,\
			'../dati/L_10/N_6/U_-1.0/bb_0.0005/fidelity_cat_s.dat' w l ls 1 t 'V_0 = 0.0005' ,\
			'../dati/L_10/N_6/U_-1.0/bb_0.0003/fidelity_cat_s.dat' w l ls 2 t 'V_0 = 0.0003' ,\
			'../dati/L_10/N_6/U_-1.0/bb_0.0001/fidelity_cat_s.dat' w l ls 3 t 'V_0 = 0.0001' ,\
			1 w l ls 9 notitle

#	set title 'N=6 L=10 U/J=-1 ANTISYMMETRIC CAT'
#	plot 	'../dati/L_10/N_6/U_-1.0/bb_0.05/fidelity_cat_a.dat' 	w l ls 2 t 'V_0 = 0.05' ,\
			'../dati/L_10/N_6/U_-1.0/bb_0.03/fidelity_cat_a.dat' 	w l ls 3 t 'V_0 = 0.03' ,\
			'../dati/L_10/N_6/U_-1.0/bb_0.01/fidelity_cat_a.dat' 	w l ls 4 t 'V_0 = 0.01' ,\
			'../dati/L_10/N_6/U_-1.0/bb_0.007/fidelity_cat_a.dat' w l ls 5 t 'V_0 = 0.007' ,\
			'../dati/L_10/N_6/U_-1.0/bb_0.005/fidelity_cat_a.dat' w l ls 6 t 'V_0 = 0.005' ,\
			'../dati/L_10/N_6/U_-1.0/bb_0.003/fidelity_cat_a.dat' w l ls 7 t 'V_0 = 0.003' ,\
			'../dati/L_10/N_6/U_-1.0/bb_0.001/fidelity_cat_a.dat' w l ls 8 t 'V_0 = 0.001' ,\
			'../dati/L_10/N_6/U_-1.0/bb_0.0007/fidelity_cat_a.dat' w l ls 9 t 'V_0 = 0.0007' ,\
			'../dati/L_10/N_6/U_-1.0/bb_0.0005/fidelity_cat_a.dat' w l ls 1 t 'V_0 = 0.0005' ,\
			'../dati/L_10/N_6/U_-1.0/bb_0.0003/fidelity_cat_a.dat' w l ls 2 t 'V_0 = 0.0003' ,\
			'../dati/L_10/N_6/U_-1.0/bb_0.0001/fidelity_cat_a.dat' w l ls 3 t 'V_0 = 0.0001' ,\
			1 w l ls 9 notitle

#	set title 'N=6 L=10 U/J=-1 CURRENT'
#	plot 	'../dati/L_10/N_6/U_-1.0/bb_0.05/corrente.dat' 	w l ls 2 t 'V_0 = 0.05' ,\
			'../dati/L_10/N_6/U_-1.0/bb_0.03/corrente.dat' 	w l ls 3 t 'V_0 = 0.03' ,\
			'../dati/L_10/N_6/U_-1.0/bb_0.01/corrente.dat' 	w l ls 4 t 'V_0 = 0.01' ,\
			'../dati/L_10/N_6/U_-1.0/bb_0.007/corrente.dat' w l ls 5 t 'V_0 = 0.007' ,\
			'../dati/L_10/N_6/U_-1.0/bb_0.005/corrente.dat' w l ls 6 t 'V_0 = 0.005' ,\
			'../dati/L_10/N_6/U_-1.0/bb_0.003/corrente.dat' w l ls 7 t 'V_0 = 0.003' ,\
			'../dati/L_10/N_6/U_-1.0/bb_0.001/corrente.dat' w l ls 8 t 'V_0 = 0.001' ,\
			'../dati/L_10/N_6/U_-1.0/bb_0.0007/corrente.dat' w l ls 9 t 'V_0 = 0.0007' ,\
			'../dati/L_10/N_6/U_-1.0/bb_0.0005/corrente.dat' w l ls 1 t 'V_0 = 0.0005' ,\
			'../dati/L_10/N_6/U_-1.0/bb_0.0003/corrente.dat' w l ls 2 t 'V_0 = 0.0003' ,\
			'../dati/L_10/N_6/U_-1.0/bb_0.0001/corrente.dat' w l ls 3 t 'V_0 = 0.0001' ,\
			1 w l ls 9 notitle







	unset multiplot
	set output	

	
	!pstopdf total.eps
	!rm total.eps



