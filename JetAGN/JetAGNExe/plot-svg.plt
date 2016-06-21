set ylabel 'Log(E_{NT} / erg s^{-1})'
set xlabel 'Log(z / pc)'
#set xrange[5.99:13]
#set yrange[7.7:18.0]
plot './'.filename.'.txt' u 1:11 w line

# plot './'.filename.'.txt' u 1:11 w line

# set terminal wxt enhanced

set terminal svg size 800,600 fname 'Verdana' fsize 10
set output filename.'.svg'

replot
set output