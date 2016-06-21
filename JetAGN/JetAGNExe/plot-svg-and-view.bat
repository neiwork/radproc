@echo off
gnuplot.exe -e "filename='%1'" ../plot-svg-all-t.plt
start chrome %~dpnx1.svg