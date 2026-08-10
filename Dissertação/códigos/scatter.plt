set title "Mode Shape"
set xlabel "X"
set ylabel "Y"
set palette defined (0 "blue", 0.5 "white", 1 "red")
set colorbox
plot 'mode1.dat' using 1:2:3 with points pt 7 ps 1.5 palette title 'u(x,y)'