load 'plt/viridispalette.pal'

set terminal pngcairo size 1000,1000 enhanced font 'Verdana,14'
set output 'plt/surfaceplot.png'

datafile='/home/thomasbjarne/Dokumenter/fortran/master_degree_project/modules/data/datafile1.txt'

set hidden3d
set view map

unset key
set size square

#set title 'q(x,t) with naive central flux and rk4'

set xlabel 'index i'
set ylabel 'time step'

set autoscale xfix
set autoscale yfix

splot datafile matrix w image