data = filename.'.txt'

set ylabel 'Log(E_{NT} / erg s^{-1})'
set xlabel 'Log(z / pc)'

#set xrange[5.99:13]
#set yrange[7.7:18.0]

# plot data u 1:11 w line

set terminal svg size 800,600 fname 'Verdana' fsize 10
set output filename.'.svg'

plot for [i=2:11] data using 1:i with line title columnhead
