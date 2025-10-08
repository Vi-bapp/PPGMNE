# Define o terminal de saída para PNG (imagem)
set terminal png size 800,600
set output 'convergence_l2_meshes.png'

# Define as escalas dos eixos como logarítmicas
set logscale x
set logscale y

# Define os rótulos dos eixos
set xlabel "Mesh size (h)"
set ylabel "L² Error"

# Define o título do gráfico
set title "L² Convergence: FEM vs. VEM - Structured and Unstructured Meshes"

# Define o grid para melhor visualização
set grid

# Define o formato dos tiques nos eixos logarítmicos
set format x "10^{%L}"
set format y "10^{%L}"

# Plota as quatro curvas em um único gráfico com cores diferentes
plot 'l2_errors_2.dat' using 1:4 with linespoints pt 7 ps 1.5 lc rgb "blue" title "FEM (Unstructured)", \
     '' using 1:2 with linespoints pt 7 ps 1.5 lc rgb "red" title "FEM (Structured)", \
     '' using 1:5 with linespoints pt 5 ps 1.5 lc rgb "green" title "VEM (Unstructured)", \
     '' using 1:3 with linespoints pt 5 ps 1.5 lc rgb "orange" title "VEM (Structured)"