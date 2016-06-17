@echo off
pgnuplot.exe -e "filename='%1'" ../plot-svg.gp
start chrome ./%1.svg