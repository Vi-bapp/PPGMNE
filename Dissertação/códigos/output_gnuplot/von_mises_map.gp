reset
set terminal pngcairo size 800,600 enhanced font 'Arial,12'
set output 'output_gnuplot/von_mises_map.png'
set title 'Tensão Equivalente de Von Mises (Pa)'
set xlabel 'X'
set ylabel 'Y'
set size ratio -1
set grid lc rgb '#e0e0e0' dt 2
set palette defined (-2 'blue', -1 'cyan', 0 'white', 1 'yellow', 2 'red')
set colorbox
plot 'output_gnuplot/von_mises_map.dat' using 1:2:3 with filledcurves lc palette title ''
