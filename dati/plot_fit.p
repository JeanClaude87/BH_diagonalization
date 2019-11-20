set autoscale


	set terminal postscript fontfile add '/Users/naldesi/Type1/sfrm1200.pfb' \
							fontfile add '/Users/naldesi/Type1/sfti1200.pfb' \
							fontfile add '/Users/naldesi/Type1/cmmi10.pfb'



!sort -n -k2  datafit_2.dat > data_fit2b.dat
!mv data_fit2b.dat datafit_2.dat 
!putblank datafit_2.dat 2 2
					
!sort -n -k2  datafit_3.dat > data_fit3b.dat
!mv data_fit3b.dat datafit_3.dat 
!putblank datafit_3.dat 2 2

!sort -n -k2  datafit_4.dat > data_fit4b.dat
!mv data_fit4b.dat datafit_4.dat 
!putblank datafit_4.dat 2 2

!sort -n -k2  datafit_5.dat > data_fit5b.dat
!mv data_fit5b.dat datafit_5.dat 
!putblank datafit_5.dat 2 2

							
	set terminal postscript	eps color enhanced font "SFRM1200" 22	size 10,5
	set border lw 2 front

	set output	'fit_dat.eps'		
 	
	set multiplot layout 2,2

	set key font "SFRM1200,22" vertical samplen 1.6 spacing 1.2

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

	set xlabel 'U'	offset 0.0,0.0	norotate 	font "SFTI1200,28" 
	set xtic 	offset 0.0,0.0	 mirror		font "SFRM1200,22"
	set mxtic			

	set ylabel 'csi'	offset 0.0,0.0 	rotate 		font "SFTI1200,28" 
	set ytic		offset 0.0,0.0	 mirror		font "SFRM1200,22"
	set mytic	

	set xrange[*:*]
	set yrange[0.005:5]

	set log 

plot 	'datafit_2.dat' index 0  u 3:4 w lp ls 2 t 'L=8',\
		''				index 1  u 3:4 w lp ls 3 t 'L=12',\
		''				index 2  u 3:4 w lp ls 4 t 'L=16',\
		''				index 3  u 3:4 w lp ls 5 t 'L=20',\
		''				index 4  u 3:4 w lp ls 6 t 'L=24',\
		''				index 5  u 3:4 w lp ls 7 t 'L=28',\
		''				index 6  u 3:4 w lp ls 8 t 'L=32',\
		''				index 7  u 3:4 w lp ls 9 t 'L=36'

plot 	'datafit_3.dat' index 0  u 3:4 w lp ls 2 t 'L=8',\
		''				index 1  u 3:4 w lp ls 3 t 'L=12',\
		''				index 2  u 3:4 w lp ls 4 t 'L=16',\
		''				index 3  u 3:4 w lp ls 5 t 'L=20',\
		''				index 4  u 3:4 w lp ls 6 t 'L=24',\
		''				index 5  u 3:4 w lp ls 7 t 'L=28',\
		''				index 6  u 3:4 w lp ls 8 t 'L=32',\
		''				index 7  u 3:4 w lp ls 9 t 'L=36'

plot 	'datafit_4.dat' index 0  u 3:4 w lp ls 2 t 'L=8',\
		''				index 1  u 3:4 w lp ls 3 t 'L=12',\
		''				index 2  u 3:4 w lp ls 4 t 'L=16',\
		''				index 3  u 3:4 w lp ls 5 t 'L=20',\
		''				index 4  u 3:4 w lp ls 6 t 'L=24'

plot 	'datafit_5.dat' index 0  u 3:4 w lp ls 2 t 'L=14'

	unset multiplot

	!pstopdf fit_dat.eps
	!rm fit_dat.eps



	set terminal postscript	eps color enhanced font "SFRM1200" 22	size 8,3
	set output	'fit_dat_tot.eps'		
	set multiplot layout 1,2

unset key

plot 	'datafit_2.dat' index 0  u 3:4 w lp ls 2 t 'L=8',\
		''				index 1  u 3:4 w lp ls 3 t 'L=12',\
		''				index 2  u 3:4 w lp ls 4 t 'L=16',\
		''				index 3  u 3:4 w lp ls 5 t 'L=20',\
		''				index 4  u 3:4 w lp ls 6 t 'L=24',\
		''				index 5  u 3:4 w lp ls 7 t 'L=28',\
		''				index 6  u 3:4 w lp ls 8 t 'L=32',\
		''				index 7  u 3:4 w lp ls 9 t 'L=36',\
	 	'datafit_3.dat' index 0  u 3:4 w lp ls 2 t 'L=8',\
		''				index 1  u 3:4 w lp ls 3 t 'L=12',\
		''				index 2  u 3:4 w lp ls 4 t 'L=16',\
		''				index 3  u 3:4 w lp ls 5 t 'L=20',\
		''				index 4  u 3:4 w lp ls 6 t 'L=24',\
		''				index 5  u 3:4 w lp ls 7 t 'L=28',\
		''				index 6  u 3:4 w lp ls 8 t 'L=32',\
		''				index 7  u 3:4 w lp ls 9 t 'L=36',\
	 	'datafit_4.dat' index 0  u 3:4 w lp ls 2 t 'L=8',\
		''				index 1  u 3:4 w lp ls 3 t 'L=12',\
		''				index 2  u 3:4 w lp ls 4 t 'L=16',\
		''				index 3  u 3:4 w lp ls 5 t 'L=20',\
		''				index 4  u 3:4 w lp ls 6 t 'L=24' 

	unset log 

	set xlabel 'U'	offset 0.0,0.0	norotate 	font "SFTI1200,28" 
	set xtic 		offset 0.0,0.0	 mirror		font "SFRM1200,22"
	set mxtic			

	set ylabel 'csi'	offset 0.0,0.0 	rotate 		font "SFTI1200,28" 
	set ytic	0.01	offset 0.0,0.0	 mirror		font "SFRM1200,22"
	set mytic	

set key
set yrange [0.085:0.115]
set xrange [0.45:0.75]


plot 	'datafit_2.dat' index 7  u 3:4 w lp ls 1 t '2,L=36',\
	 	'datafit_3.dat' index 5  u 3:4 w lp ls 2 t '3,L=28',\
	 	'datafit_4.dat' index 3  u 3:4 w lp ls 3 t '4,L=20',\
	 	'datafit_5.dat' index 1  u 3:4 w lp ls 4 t '5,L=14',\
	 	'datafit_2.dat' index 6  u 3:4 w lp ls 5 t '2,L=32',\
	 	'datafit_3.dat' index 4  u 3:4:5 w e ls 6 t '3,L=24',\
	 	'datafit_4.dat' index 2  u 3:4:5 w e ls 7 t '4,L=16',\
	 	'datafit_5.dat' index 0  u 3:4:5 w e ls 4 t '5,L=10',\
	 	0.10

	!pstopdf fit_dat_tot.eps
	!rm fit_dat_tot.eps





