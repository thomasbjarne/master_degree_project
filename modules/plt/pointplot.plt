datafile='/home/thomasbjarne/Dokumenter/fortran/master_degree_project/modules/data/datafile1.txt'

set terminal pngcairo size 1280,720 enhanced font 'Verdana,12'
set output 'plt/pointplot.png'

set size square

unset key
#unset border
#unset tics


plot datafile w lp