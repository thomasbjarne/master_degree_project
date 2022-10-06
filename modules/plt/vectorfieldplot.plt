set terminal pngcairo size 1000,1000 enhanced font 'Verdana,14'
set output 'plt/vectorfieldplot.png'

datafile='/home/thomasbjarne/Dokumenter/fortran/master_degree_project/modules/data/datafile1.txt'

set hidden3d
set view map

unset key
set size square

#set title 'q(x,t) with naive central flux and rk4'

set xlabel 'x'
set ylabel 'y'

set autoscale xfix
set autoscale yfix

plot datafile u 1:2:3:4 w vectors head filled lt 2