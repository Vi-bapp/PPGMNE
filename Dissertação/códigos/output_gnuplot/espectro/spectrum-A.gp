# =====================================================================
# ESPECTRO DE FREQUENCIAS - SCRIPT GENERICO COM CABECALHO NO CSV
# CSV: campos separados por ';', ponto decimal.
# Coluna 1 = referencia; colunas 2 a 7 = series numericas.
# As legendas sao lidas automaticamente do cabecalho do CSV.
# =====================================================================

# ------------------------- CONFIGURACAO -------------------------------
FICHEIRO = "dados_excel_V.csv"
SAIDA    = "espectro_frequencias_Vb.png"

COL_REF = 1

# Escolha as colunas que deseja plotar (indices das colunas no CSV).
# Primeira ordem:
COLUNAS = "2 3 4"

# Segunda ordem:
# COLUNAS = "5 6 7"

# Comparar ordens com 4 elementos:
# COLUNAS = "2 5"

# Limites dos eixos
XMIN = 1
XMAX = 10
YMIN = 0
YMAX_AUTO = 1       # 1 = limite superior automatico; 0 = usar YMAX
YMAX = 3.0

TITULO_GRAFICO = "Espectro de Frequencias (Razao em Relacao a Referencia)"
ROTULO_X = "Modo (n)"
ROTULO_Y = "Razao (ω_num / ω_ref)"
LARGURA = 1000
ALTURA = 600
# ----------------------------------------------------------------------

set datafile separator ";"
set datafile columnheaders

set terminal pngcairo size LARGURA,ALTURA font "Arial,12"
set output SAIDA
set title TITULO_GRAFICO
set xlabel ROTULO_X
set ylabel ROTULO_Y
set grid xtics ytics
set key top left box horizontal
set xrange [XMIN:XMAX]
if (YMAX_AUTO) {
    set yrange [YMIN:*]
} else {
    set yrange [YMIN:YMAX]
}
set xtics 1

N = words(COLUNAS)

# column(0) fornece a posicao do registro. Com este CSV, o primeiro
# registro numerico deve ser plotado no modo 1 (sem somar +1).
plot for [j=1:N] FICHEIRO \
    using (column(0)+1):(column(int(word(COLUNAS,j)))/column(COL_REF)) \
    with linespoints pointtype (j+6) pointsize 0.8 linewidth 2 \
    title columnhead(int(word(COLUNAS,j)) + 0)

unset output
