#!/usr/bin/env gnuplot

# set the output as png
set term pdfcairo

# Input file contains comma-separated values fields
set output "output_data/Genergy.pdf"
set title "Gradient of energy"
set xlabel "variational parameter b"
set ylabel "Gradeint of energy"
set grid
# set xrange[:0.08]
plot "data/Genergy.dat" using 1:2:3 with yerrorbars pointtype 7 pointsize 0.3
set output