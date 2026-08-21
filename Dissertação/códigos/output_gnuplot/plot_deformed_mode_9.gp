reset
set terminal pngcairo size 800,600 enhanced font 'Arial,12'
set output 'output_gnuplot/deformed_mode_9.png'
set title 'Malha Deformada - Modo 9'
set xlabel 'X'
set ylabel 'Y'
set size ratio -1
set grid lc rgb '#e0e0e0' dt 2
set key top right
plot 'output_gnuplot/mesh_original.gp.dat' with lines lw 1.2 lc rgb '#c0c0c0' dt 2 title 'Original', 'output_gnuplot/mode_9_mesh_def.gp.dat' with lines lw 2.0 lc rgb '#0060ad' title 'Deformada'
