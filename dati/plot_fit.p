set autoscale


	set terminal postscript fontfile add '/Users/naldesi/Type1/sfrm1200.pfb' \
							fontfile add '/Users/naldesi/Type1/sfti1200.pfb' \
							fontfile add '/Users/naldesi/Type1/cmmi10.pfb'
							
							
	set terminal postscript	eps color enhanced font "SFRM1200" 22	size 10,3
	set border lw 2 front

	set output	'fit_dat.eps'		
 	
	set multiplot layout 1,3

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
	set yrange[*:*]

	set log 

plot 	'datafit_2.dat' index 0  u ($3*$2):4 w lp ls 2 t 'L=10',\
		''				index 1  u ($3*$2):4 w lp ls 3 t 'L=14',\
		''				index 2  u ($3*$2):4 w lp ls 4 t 'L=18',\
		''				index 3  u ($3*$2):4 w lp ls 5 t 'L=22',\
		''				index 4  u ($3*$2):4 w lp ls 6 t 'L=26',\
		''				index 5  u ($3*$2):4 w lp ls 7 t 'L=30',\
		''				index 6  u ($3*$2):4 w lp ls 8 t 'L=34',\
		''				index 7  u ($3*$2):4 w lp ls 9 t 'L=38'

plot 	'datafit_3.dat' index 0  u ($3*$2):4 w lp ls 2 t 'L=10',\
		''				index 1  u ($3*$2):4 w lp ls 3 t 'L=14',\
		''				index 2  u ($3*$2):4 w lp ls 4 t 'L=18',\
		''				index 3  u ($3*$2):4 w lp ls 5 t 'L=22',\
		''				index 0  u ($3*$2):4:5 w e ls 2 notitle,\
		''				index 1  u ($3*$2):4:5 w e ls 3 notitle,\
		''				index 2  u ($3*$2):4:5 w e ls 4 notitle,\
		''				index 3  u ($3*$2):4:5 w e ls 5 notitle

plot 	'datafit_4.dat' index 3  u ($3*$2):4 w lp ls 2 t 'L=8',\
		''				index 0  u ($3*$2):4 w lp ls 3 t 'L=12',\
		''				index 1  u ($3*$2):4 w lp ls 4 t 'L=16',\
		''				index 2  u ($3*$2):4 w lp ls 5 t 'L=20',\
		''				index 3  u ($3*$2):4:5 w e ls 2 notitle,\
		''				index 0  u ($3*$2):4:5 w e ls 3 notitle,\
		''				index 1  u ($3*$2):4:5 w e ls 4 notitle,\
		''				index 2  u ($3*$2):4:5 w e ls 5 notitle









	!pstopdf fit_dat.eps
	!rm fit_dat.eps

