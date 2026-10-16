reset
set terminal pngcairo size 900,500 enhanced font 'Arial,12'
set output 'output_gnuplot/espectro/vem2-2_locking.png'
set title 'Indices de locking (VEM) vs ordem modal adimensional'
set xlabel 'Ordem modal adimensional i/N'
set ylabel 'r'
set yrange [0:1]
set grid lc rgb '#e0e0e0' dt 2
set key top right
plot 'output_gnuplot/espectro/vem2-2_locking.dat' using 1:2 with linespoints lw 2 pt 7 lc rgb '#d9534f' title 'r_{vol}', \
     'output_gnuplot/espectro/vem2-2_locking.dat' using 1:3 with linespoints lw 2 pt 7 lc rgb '#0060ad' title 'r_{shear}'
