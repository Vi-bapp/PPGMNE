# ===============================================================
# PLOT DEFORMADA - MEV PÓRTICOS PLANOS 2D
# Arquivo: plot_deformada.gp
# Execute: gnuplot plot_deformada.gp
# COMPATÍVEL COM TODAS AS VERSÕES
# ===============================================================

# --- Saída PNG ---
set terminal pngcairo enhanced font 'Arial,12' size 900,600
set output 'deformada_mev.png'

set title "MEV - Malha Deformada (Pórticos Planos 2D)" font "Arial,14"
set xlabel "x (m)" font "Arial,12"
set ylabel "y (m)" font "Arial,12"
set grid
set key top left box

# --- Escala de deslocamentos ---
scale = 50.0

# --- Legenda personalizada ---
legend_text = sprintf("Deformada (×%.0f)", scale)

# --- Plot 1: Malha original (preta) ---
plot 'solution_vem_plane_frame.dat' using 2:3 with lines lc rgb "black" lw 2 title "Malha Original", \
     '' using 2:3 with points pt 7 ps 1.5 lc rgb "black" notitle

# --- Plot 2: Malha deformada (azul) ---
replot 'solution_vem_plane_frame.dat' using 7:8 with lines lc rgb "blue" lw 2.5 title legend_text, \
       '' using 7:8 with points pt 7 ps 1.8 lc rgb "blue" notitle

# --- Salva em PDF ---
set terminal pdfcairo enhanced font 'Arial,12'
set output 'deformada_mev.pdf'
replot

# --- Finaliza ---
set output