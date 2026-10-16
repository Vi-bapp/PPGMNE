reset
set terminal pngcairo size 1000,550 enhanced font 'Arial,12'
set output 'output_gnuplot/espectro/vem1-3_locking.png'
set title 'Indices de locking (VEM) vs ordem modal adimensional'
set xlabel 'Ordem modal adimensional i/N'
set ylabel 'Indice de locking r'
set yrange [0:1.05]
set grid lc rgb '#e0e0e0' dt 2
set key top right
set border linewidth 1.5
plot 'output_gnuplot/espectro/vem1-3_locking.dat' using 2:3 with lines lw 2 lc rgb '#d9534f' title 'r_{vol}', \
     'output_gnuplot/espectro/vem1-3_locking.dat' using 2:4 with lines lw 2 lc rgb '#0060ad' title 'r_{shear}'
