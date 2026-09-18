reset
set terminal pngcairo size 900,700 enhanced font 'Arial,12'
set output 'output_gnuplot/sigma_xx_map.png'
set title 'Tensao Normal XX (Pa)'
set xlabel 'X'
set ylabel 'Y'
set size ratio -1
set grid lc rgb '#e0e0e0' dt 2
set palette defined (0 'blue', 0.25 'cyan', 0.5 'green', 0.75 'yellow', 1 'red')
set colorbox
set cbrange [ -2.725431E+01:  2.963846E+01]
set cblabel 'field'
plot 'output_gnuplot/sigma_xx_map.dat' using 1:2:3 with filledcurves closed lc palette title ''
