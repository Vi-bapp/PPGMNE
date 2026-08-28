reset
set terminal pngcairo size 800,400 enhanced font 'Arial,12'
set output 'output_gnuplot/time_disp_mode_3_dof_10.png'
set title 'Resposta no Tempo - Modo 3 (GL 10)'
set xlabel 'Tempo (s)'
set ylabel 'Deslocamento u_{10}(t)'
set grid lc rgb '#e0e0e0' dt 2
unset key
plot 'output_gnuplot/time_disp_mode_3_dof_10.dat' using 1:2 with lines lw 2 lc rgb '#d9534f'
