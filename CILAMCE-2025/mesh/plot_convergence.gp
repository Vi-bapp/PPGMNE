# Define o terminal de saída para PNG (imagem)
set terminal pngcairo enhanced font "Arial,12"
set output 'convergence_l2.png'

# Define as escalas dos eixos como logarítmicas
set logscale x
set logscale y

# Define os rótulos dos eixos
set xlabel "Mesh size (h)"
set ylabel "L² Error"

# Define o título do gráfico
set title "L² convergence: FEM vs. VEM"

# Define o grid para melhor visualização
set grid

# Define o formato dos tiques nos eixos logarítmicos
set format x "10^{%L}"
set format y "10^{%L}"

# Define as funções de referência para convergência de primeira ordem (slope = 1)
# Ajuste as constantes 'C_fem' e 'C_vem' para que as linhas teóricas
# se ajustem visualmente aos seus primeiros pontos de dados.
# Para FEM (0.25, 0.054936): log(0.054936) = log(C_fem) + log(0.25) => -2.90 = log(C_fem) - 1.386 => log(C_fem) = -1.514 => C_fem = exp(-1.514) = 0.22
# Para VEM (0.25, 0.079848): log(0.079848) = log(C_vem) + log(0.25) => -2.527 = log(C_vem) - 1.386 => log(C_vem) = -1.141 => C_vem = exp(-1.141) = 0.319
C_fem = 0.22
C_vem = 0.32
f_fem(x) = C_fem * x**1 # Convergência de primeira ordem (FEM)
f_vem(x) = C_vem * x**1 # Convergência de primeira ordem (VEM)

# Plota os dados e as linhas de referência
plot 'l2_errors.dat' using 1:2 with linespoints pt 7 ps 1.5 lc rgb "blue" title "(FEM)", \
     '' using 1:3 with linespoints pt 5 ps 1.5 lc rgb "green" title "(VEM)", \

# Define os estilos das linhas teóricas
set style line 1 lt 1 lw 2 # Linha sólida, espessura 2
set style line 2 lt 2 lw 2 # Linha tracejada, espessura 2