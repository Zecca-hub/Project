#!/usr/bin/env gnuplot

# set the output as png
set term pdfcairo
# Input file contains comma-separated values fields
set output "output_data/Eac.pdf"
set title "Energy autocorrelation"
set xlabel "Step"
set ylabel "Eac"
# set log y
# set xrange [0:550]
# set yrange [-100:400]
l0=1
l1=3.05462e-06
l2=0
set grid
plot "data/Eac.dat" using 1:2 with lines notitle
# , l0*exp(-l1*x)+l2 with lines notitle
# plot "output_data/Qautocorr.dat" using 2:3 with lines notitle
# write to file
set output