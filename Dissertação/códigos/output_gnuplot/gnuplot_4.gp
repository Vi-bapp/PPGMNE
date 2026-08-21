set terminal pngcairo size 800,600 enhanced font 'Arial,12'
set output 'grafico_4_eficiencia.png'

set logscale y
set format y "10^{%L}"
set xlabel "Tempo de CPU total (s)" font ",14" bold
set ylabel "Erro Relativo da 1ª Freq. (%)" font ",14" bold
set grid lc rgb '#e0e0e0' dt 2
set key top right

plot "vem_dados_globais.dat" using 4:2 with points lc rgb '#0060ad' pt 7 ps 2 title "VEM", \
     "mefg_tempo_erro.dat" using 1:2 with points lc rgb '#dd181f' pt 5 ps 2 title "MEFG"