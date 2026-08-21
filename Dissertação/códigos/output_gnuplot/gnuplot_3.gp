set terminal pngcairo size 800,600 enhanced font 'Arial,12'
set output 'grafico_3_espectro.png'

set xlabel "Ordem Normalizada (i/N)" font ",14" bold
set ylabel "Frequência Normalizada (w / w_{ref})" font ",14" bold
set xrange [0:1]
set yrange [0.8:2.0]
set grid lc rgb '#e0e0e0' dt 2

set arrow from 0,1 to 1,1 nohead lc rgb "black" dt 2 lw 2

plot "vem_espectro.dat" using 1:2 with lines lc rgb '#0060ad' lw 3 title "VEM", \
     "mefg_espectro.dat" using 1:2 with lines lc rgb '#dd181f' lw 3 title "MEFG-2 (Custódio)"