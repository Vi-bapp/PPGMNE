reset
set terminal pngcairo size 800,600 enhanced font 'Arial,12'
set output 'modo_2.png'
set title 'Superfície do Modo 2'
set xlabel 'X'
set ylabel 'Y'
set zlabel 'Deslocamento'
set dgrid3d 50,50 qnorm 2
set pm3d
set hidden3d
set view 60, 30
set palette defined (0 'blue', 1 'cyan', 2 'green', 3 'yellow', 4 'red')
splot 'modo_2.dat' with pm3d notitle
