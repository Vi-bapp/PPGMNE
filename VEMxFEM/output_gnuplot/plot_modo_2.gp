reset
set terminal pngcairo size 800,600 enhanced font 'Arial,12'
set output 'output_gnuplot/modo_2.png'
set title 'Campo de Deslocamento |u| - Modo 2'
set xlabel 'X'
set ylabel 'Y'
set size ratio -1
set dgrid3d 60,60 qnorm 2
set pm3d map
set palette defined (0 'blue', 1 'cyan', 2 'green', 3 'yellow', 4 'red')
splot 'output_gnuplot/modo_2.dat' using 1:2:3 with pm3d notitle
