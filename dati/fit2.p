set autoscale
set terminal postscript fontfile add '/Users/naldesi/Type1/sfrm1200.pfb'\
						fontfile add '/Users/naldesi/Type1/sfti1200.pfb' 			
set terminal postscript	eps color enhanced font "SFRM1200" 32	size 5,3


first(x) = ($0 > 0 ? base : base = x)

lista = system('ls N_5/corr*.dat')

set fit errorvariables
set print 'datafit_5.dat'

f(x) = a*( exp(-x/b) )

do for [xx in lista] { #}

	dati = system(sprintf('echo %s | sed -e "s/\.dat/ /g" -e "s/N/ /g" -e "s/corr/ /g" -e "s/L/ /g" | tr "_" " " | tr "/" " " | tr "\-U" " " | sed -e "s/  / /g" ', xx))
	n0 = system(sprintf("echo %s | awk '{print $1}' ", dati))
	l0 = system(sprintf("echo %s | awk '{print $2}' ", dati))
	u0 = system(sprintf("echo %s | awk '{print $3}' ", dati))

 	fit [0.001:0.4] f(x) xx u ($2/$1):3 via a, b  

	print n0, ' ', l0, ' ', u0, '  ', b, '  ', b_err
	
	outfile = system(sprintf('echo %s | sed -e "s/\.dat/\.eps/g" ', xx))


set output outfile #'ciao.eps'

	set log y

	set xrange[*:0.5]
	set yrange[*:*]

 	plot 	xx u ($2/$1):3 w lp ls 1, f(x) ls 2 
		
 set output



	}

quit

