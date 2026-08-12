set terminal pngcairo size 1200,800 font "Arial,12"
set grid
set key top right
set size ratio -1

# Configuração do mapa de cores para magnitude do deslocamento
set palette defined (0 "blue", 1 "cyan", 2 "green", 3 "yellow", 4 "red")
set colorbox
set cblabel "Displacement Magnitude |u|"

# 1. Plot da Malha Original vs Malha Deformada (Modo 1)
set output "mode_1_deformation.png"
set title "VEM Vibration Mode 1 - Deformed Mesh vs Original Mesh"
plot "mesh_original.gp.dat" with lines dt 2 lc rgb "gray50" title "Original Mesh", \
     "mode_1_mesh_def.gp.dat" with lines lw 2 lc rgb "navy" title "Deformed Mesh"

# 2. Plot de Campo de Vetores e Contorno de Deslocamento
set output "mode_1_field.png"
set title "VEM Vibration Mode 1 - Displacement Vector Field and Magnitude"
plot "mode_1_mesh_def.gp.dat" with lines lw 1.5 lc rgb "black" notitle, \
     "mode_1_data.gp.dat" using 2:3:4:5:6 with vectors head filled lt 2 palette lw 1.5 title "Displacement Vectors", \
     "mode_1_data.gp.dat" using 7:8:6 with points pt 7 ps 1.5 palette title "Node Magnitude"