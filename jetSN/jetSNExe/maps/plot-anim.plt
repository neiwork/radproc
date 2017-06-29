data = 'Mbh.txt'

set key autotitle columnhead # si no stats falla
stats data nooutput
COLUMNS = int(STATS_columns)

set ylabel 'Log(M / Ms)'
set xlabel 'Log(N)'

set terminal gif animate delay 10
set output 'Nmt.gif'

do for [i=2:COLUMNS] {
	plot data using 1:i with line title columnhead
}