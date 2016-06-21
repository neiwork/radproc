data = filename.'.txt'

set ylabel 'Log(E_{NT} / erg s^{-1})'
set xlabel 'Log(z / pc)'

set terminal png 
# size 800,600 fname 'Verdana' fsize 10

do for [i=2:11] {
	outputfile = sprintf(filename.'-anim-%02d.png',i-1)
	set output outputfile
	plot data using 1:i with line
}