#!/usr/bin/gnuplot
set termoption dash
set terminal pdf enhanced dashed color font 'Helvetica-bold,20'
set datafile missing "?"

set output 'QE.pdf'

set xrange [-1:7]    
set xtics rotate by 45 right
set ylabel "QuasiEntropy" #offset 2,0
set logscale y 10

set for [i=1:5] linetype i dt i

plot 'QE.data' using ($2):xticlabels(1) w linespoints ps 3 pt 7 dt 1 lw 4 lc rgb 'black' notitle ,\
