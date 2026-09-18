reset
set terminal pngcairo size 900,700 enhanced font 'Arial,12'
set output 'output_gnuplot/tau_xy_map.png'
set title 'Tensao de Cisalhamento XY (Pa)'
set xlabel 'X'
set ylabel 'Y'
set size ratio -1
set grid lc rgb '#e0e0e0' dt 2
set palette defined (0 'blue', 0.25 'cyan', 0.5 'green', 0.75 'yellow', 1 'red')
set colorbox
set cbrange [ -2.879372E+00:  1.521659E+01]
set cblabel 'field'
plot 'output_gnuplot/tau_xy_map.dat' using 1:2:3 with filledcurves closed lc palette title ''
