reset
set terminal gif animate delay 8 loop 0 size 800,600 enhanced font 'Arial,10'
set output 'output_gnuplot/animacao_vem.gif'
set title 'Evolucao Temporal da Malha (VEM)'
set xlabel 'X'
set ylabel 'Y'
set xrange [  5.639205E+03:  5.642005E+03]
set yrange [  1.246110E+03:  1.248910E+03]
set size ratio -1
set lmargin at screen 0.12
set rmargin at screen 0.92
set bmargin at screen 0.12
set tmargin at screen 0.88
set grid lc rgb '#e0e0e0' dt 2
unset key
do for [i=1:40] {
    frame_file = sprintf('output_gnuplot/deformed_frame_%04d.gp.dat', i)
    plot 'output_gnuplot/mesh_original.gp.dat' with lines lw 1 lc rgb '#a0a0a0' dt 2, \
         frame_file with lines lw 2 lc rgb '#0060ad'
}
set output
