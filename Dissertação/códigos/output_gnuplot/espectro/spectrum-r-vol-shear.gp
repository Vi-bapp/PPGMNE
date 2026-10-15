# =====================================================================
# ESPECTRO DOS INDICADORES r_vol e r_shear
# Eixo x: modos normalizados no intervalo [0,1]
# CSV: separador ';', cabeçalho na primeira linha.
# Colunas esperadas:
#   1 = Modo | 2 = omega_rad_s | 3 = f_Hz
#   4 = r_vol | 5 = r_shear | 6 = lambda
# =====================================================================

# ------------------------- CONFIGURACAO -------------------------------
# Informe os arquivos CSV a comparar, separados por espacos.
ARQUIVOS = "vem1-1.csv vem1-2.csv vem1-3.csv"

# Legendas correspondentes, na mesma ordem dos arquivos.
TITULOS = "Malha_1"

SAIDA = "espectro_r_vol_r_shear.png"

TITULO_GRAFICO = "Indicadores r_vol e r_shear por modo normalizado"
ROTULO_X = "Modo normalizado, (n-1)/(N-1)"
ROTULO_Y = "Indicador"

LARGURA = 1100
ALTURA = 650

# Limites verticais: use YMAX_AUTO=1 para limite automatico.
YMIN = 0.0
YMAX_AUTO = 1
YMAX = 1.0
# ----------------------------------------------------------------------

set datafile separator ";"
set datafile columnheaders

set terminal pngcairo size LARGURA,ALTURA font "Arial,12"
set output SAIDA
set title TITULO_GRAFICO
set xlabel ROTULO_X
set ylabel ROTULO_Y
set xrange [0:1]
if (YMAX_AUTO) {
    set yrange [YMIN:*]
} else {
    set yrange [YMIN:YMAX]
}
set xtics 0,0.1,1
set grid xtics ytics
set key outside bottom center horizontal

NARQ = words(ARQUIVOS)

# Verifica e armazena o primeiro e ultimo modo de cada arquivo.
# A normalizacao (n-n_min)/(n_max-n_min) garante x=0 no primeiro
# modo e x=1 no ultimo modo, mesmo que a numeracao nao comece em 1.
MODOS_MIN = ""
MODOS_MAX = ""

do for [i=1:NARQ] {
    stats word(ARQUIVOS,i) using 1 nooutput
    MODOS_MIN = sprintf("%s %.12g", MODOS_MIN, STATS_min)
    MODOS_MAX = sprintf("%s %.12g", MODOS_MAX, STATS_max)
}

# Cada arquivo gera duas curvas: r_vol (coluna 4) e r_shear (coluna 5).
plot for [i=1:NARQ] word(ARQUIVOS,i) \
    using (($1-real(word(MODOS_MIN,i)))/(real(word(MODOS_MAX,i))-real(word(MODOS_MIN,i)))):4 \
    with linespoints linewidth 2 pointsize 0.7 pointtype (2*i+5) \
    title sprintf("%s - r_vol", word(TITULOS,i)), \
     for [i=1:NARQ] word(ARQUIVOS,i) \
    using (($1-real(word(MODOS_MIN,i)))/(real(word(MODOS_MAX,i))-real(word(MODOS_MIN,i)))):5 \
    with linespoints linewidth 2 pointsize 0.7 pointtype (2*i+6) dashtype 2 \
    title sprintf("%s - r_shear", word(TITULOS,i))

unset output
