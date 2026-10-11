# ==============================================================================
# SCRIPT 2: OPÇÃO B - MÚLTIPLOS ARQUIVOS DO SOLVER (w_n / w_1 vs. n / N_total)
# ==============================================================================
# Liste os nomes dos arquivos CSV gerados pelo solver e as legendas do gráfico:
ARQUIVOS = "cu4n9k1.csv cu16n25k1.csv cu64n81k1.csv cu256n289k1.csv"
TITULOS  = "cu4 cu16 cu64 cu256"

set datafile separator ";"
set decimalsign ','      # Descomente se os CSVs do solver usarem vírgula decimal

# ------------------------------------------------------------------------------
# PRÉ-CÁLCULO AUTOMÁTICO DE N_TOTAL E w_1 PARA CADA ARQUIVO DO SOLVER
# ------------------------------------------------------------------------------
num_arqs     = words(ARQUIVOS)
N_TOTAL_AUTO = ""
W1_AUTO      = ""

do for [i=1:num_arqs] {
    set xrange [*:*]
    set yrange [*:*]
    
    # Obtém total de modos (N_total) e a 1ª frequência (w_1) da malha
    stats word(ARQUIVOS, i) every ::1 using 1:2 nooutput
    
    N_TOTAL_AUTO = sprintf("%s %d", N_TOTAL_AUTO, STATS_records)
    W1_AUTO      = sprintf("%s %f", W1_AUTO, STATS_min_y)
}

# ------------------------------------------------------------------------------
# CONFIGURAÇÃO DE TERMINAL E EIXOS
# ------------------------------------------------------------------------------
set terminal pngcairo size 1000,600 font "Arial,12"
set output 'espectro_opcao_B.png'

set grid xtics ytics mxtics mytics
set key top left box horizontal

set title "Espectro Adimensional de Frequências ({/Symbol w}_n / {/Symbol w}_1)"
set xlabel "n / N_{total}"
set ylabel "{/Symbol w}_n / {/Symbol w}_1"

set xrange [0.0 : 1.0]
set yrange [0 : *]
set xtics auto

# ------------------------------------------------------------------------------
# PLOTAGEM DINÂMICA DE TODOS OS ARQUIVOS
# ------------------------------------------------------------------------------
# Coluna 1: Índice do Modo (n)
# Coluna 2: Frequência (w_n)
plot for [i=1:num_arqs] \
    word(ARQUIVOS, i) every ::1 \
    using ($1 / (word(N_TOTAL_AUTO, i) + 0.0)):($2 / (word(W1_AUTO, i) + 0.0)) \
    w l lw 2 title word(TITULOS, i)

unset output