set terminal pngcairo size 800,600 enhanced font 'Arial,12'
set output 'grafico_2_condicionamento.png'

set logscale y
set format y "10^{%L}"
set xlabel "Graus de Liberdade (G.L.)" font ",14" bold
set ylabel "Número de Condição da Matriz" font ",14" bold
set grid xtics ytics mxtics mytics lc rgb '#e0e0e0' dt 2
set key left top

set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 ps 1.5
set style line 2 lc rgb '#dd181f' lt 2 lw 2 pt 9 ps 1.5

plot "vem_dados_globais.dat" using 1:3 with linespoints ls 1 title "VEM", \
     "mefg_condicao.dat" using 1:2 with linespoints ls 2 title "MEFG"