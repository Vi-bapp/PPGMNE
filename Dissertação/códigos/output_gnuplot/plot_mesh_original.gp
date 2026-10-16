reset
set terminal pngcairo size 800,600 enhanced font 'Arial,12'
set output 'output_gnuplot/mesh_original.png'
set title 'Malha Original'
set xlabel 'X'
set ylabel 'Y'
set size ratio -1
set grid lc rgb '#e0e0e0' dt 2
unset key
plot 'output_gnuplot/mesh_original.gp.dat' with lines lw 1.5 lc rgb '#333333'
