# Define o estilo do gráfico
set style data lines
set xrange [-10:110] # Ajuste para as dimensões da sua estrutura
set yrange [-10:10]

# Loop para plotar todos os quadros
do for [i=10:1000:10] {
    filename = sprintf("dyn_step_%05d.dat", i)
    plot filename title sprintf("Tempo: Passo %d", i)
    pause 0.05
}