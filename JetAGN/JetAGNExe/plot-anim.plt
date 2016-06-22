data = filename.'.txt'

set key autotitle columnhead # si no stats falla
stats data nooutput
COLUMNS = int(STATS_columns)

set ylabel 'Log(E_{NT} / erg s^{-1})'
set xlabel 'Log(z / pc)'

set terminal gif animate delay 10
set output filename.'.gif'

do for [i=2:COLUMNS] {
	plot data using 1:i with line title columnhead
}