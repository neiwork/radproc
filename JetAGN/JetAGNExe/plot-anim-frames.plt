data = filename.'.txt'

set key autotitle columnhead # si no stats falla
stats data nooutput
COLUMNS = int(STATS_columns)

set ylabel 'Log(E_{NT} / erg s^{-1})'
set xlabel 'Log(z / pc)'

set terminal png 

do for [i=2:COLUMNS] {
	outputfile = sprintf(filename.'-anim-%02d.png',i-1)
	set output outputfile
	plot data using 1:i with line
}