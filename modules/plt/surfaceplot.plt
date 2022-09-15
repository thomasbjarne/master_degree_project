load 'plt/viridis.pal'

set terminal pngcairo size 1280,720 enhanced font 'Verdana,12'
set output 'plt/surfaceplot.png'

datafile='/home/thomasbjarne/fortran/master_thesis/modules/data/datafile1.txt'

set hidden3d
set view map

unset key
set size square

set title 'u(x,t) with rk4 and central diff'

set xlabel 'x'
set ylabel 't'

set xtics 0, 0.25, 1
set ytics 0, 0.25, 1

splot datafile nonuniform matrix u 1:2:3 w pm3d