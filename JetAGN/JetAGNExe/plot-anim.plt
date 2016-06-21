data = filename.'.txt'

if(!exist("columns")) columns=1

set ylabel 'Log(E_{NT} / erg s^{-1})'
set xlabel 'Log(z / pc)'

set terminal gif animate delay 10
set output filename.'.gif'

do for [i=2:int(columns+1)] {
	plot data using 1:i with line title columnhead
}