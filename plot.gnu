#!/usr/bin/gnuplot 

set term pdf monochrome size 15.0cm,9.0cm

set output "pressure.pdf"

set grid

set xtics font "Times-Roman, 13"
set ytics font "Times-Roman, 13"

set xlabel "x position" center
set ylabel "Pressure [-]" center
set title "Analytical vs Exact - Oblique shock reflection"

set offset graph 0.10,0.10,0.10,0.10

set border lw 2

set pointsize 0.3

set key font ",12"

plot "exact.sol" u 1:2 with linespoints lt -1 lw 0.1 pt 4 title "Analytical",\
     "mid_pressure.dat" u 1:2 with lines lw 1.5 title "Numerical"
