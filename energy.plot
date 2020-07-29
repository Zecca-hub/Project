#!/usr/bin/env gnuplot

# set the output as png
set term pdfcairo

# Input file contains comma-separated values fields
set output "output_data/Energy.pdf"
set title "Energy"
set xlabel "variational parameter b"
set ylabel "Energy"
set grid
# set xrange[:0.08]
plot "data/Energy.dat" using 1:2:3 with yerrorbars pointtype 7 pointsize 0.3
set output