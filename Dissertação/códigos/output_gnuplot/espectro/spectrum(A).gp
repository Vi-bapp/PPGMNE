# ==============================================================================
# SCRIPT 1: OPÇÃO A - RAZÃO EM RELAÇÃO À REFERÊNCIA (w_num / w_ref)
# ==============================================================================
FICHEIRO = "dados_excel.csv"
TITULOS  = "'1ª Ord. (4 el)' '1ª Ord. (16 el)' '1ª Ord. (64 el)' '2ª Ord. (4 el)' '2ª Ord. (16 el)' '2ª Ord. (64 el)'"

set datafile separator ";"
set decimalsign ','      # Descomente/Mantenha se o CSV usa vírgula decimal

# ------------------------------------------------------------------------------
# CONFIGURAÇÃO DO TERMINAL E EIXOS
# ------------------------------------------------------------------------------
set terminal pngcairo size 1000,600 font "Arial,12"
set output 'espectro_opcao_A.png'

set grid xtics ytics mxtics mytics
set key top left box horizontal

set title "Espectro de Frequências (Razão em Relação à Referência)"
set xlabel "Modo (n)"
set ylabel "Razão ({/Symbol w}_{num} / {/Symbol w}_{ref})"

set xrange [1 : 10]
set yrange [0.8 : 5.0]
set xtics 1,1,10

# ------------------------------------------------------------------------------
# PLOTAGEM DINÂMICA
# ------------------------------------------------------------------------------
# Divide a coluna da malha (3 a 8) pela coluna de referência (Coluna 2)
plot for [i=3:8] FICHEIRO every ::1 \
    using 1:(column(i) / $2) \
    w lp pt 7 ps 0.8 lw 2 title word(TITULOS, i-2)

unset output