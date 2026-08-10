set title "Mode Shape Surface"
set dgrid3d 60,60 qnorm 2    # Interpolates scattered (x,y,u) onto a 60x60 grid
set pm3d
set hidden3d
splot 'mode1.dat' using 1:2:3 with pm3d