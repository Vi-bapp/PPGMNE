reset
set terminal pngcairo size 800,400 enhanced font 'Arial,12'
set output 'output_gnuplot/time_response_mode_1.png'
set title 'Coordenada Modal q_1(t) - Modo 1'
set xlabel 'Tempo (s)'
set ylabel 'Amplitude q_n(t)'
set grid lc rgb '#e0e0e0' dt 2
unset key
plot 'output_gnuplot/time_response_mode_1.dat' using 1:2 with lines lw 2 lc rgb '#0060ad'
