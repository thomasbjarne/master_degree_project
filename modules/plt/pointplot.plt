datafile='/home/thomasbjarne/fortran/master_thesis/modules/data/datafile1.txt'

set terminal pngcairo size 1280,720 enhanced font 'Verdana,12'
set output 'plt/pointplot.png'

set size square

unset key
#unset border
#unset tics


plot datafile w lp