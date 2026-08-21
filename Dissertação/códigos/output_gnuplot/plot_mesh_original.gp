reset
set terminal pngcairo size 800,600 enhanced font 'Arial,12'
set output 'output_gnuplot/mesh_original.png'
set title 'Malha Poligonal Original (VEM)'
set xlabel 'X'
set ylabel 'Y'
set size ratio -1
set grid lc rgb '#e0e0e0' dt 2
plot 'output_gnuplot/mesh_original.gp.dat' with lines lw 1.5 lc rgb '#0060ad' notitle
