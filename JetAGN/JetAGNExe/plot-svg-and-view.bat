@echo off

if "%2"=="" (
	set METHOD=plot-svg-last-t
	set EXT=svg
) else (
	set METHOD=%2
	set EXT=%3
)

pushd ..
gnuplot.exe -e "filename='%1'" ../%METHOD%.plt
if EXIST %~dp1plots\%~n1.%EXT% rm %~dp1plots\%~n1.%EXT%
mv %~n1.%EXT% plots
echo opening %~dp1plots\%~n1.%EXT%
start chrome %~dp1plots\%~n1.%EXT%
popd