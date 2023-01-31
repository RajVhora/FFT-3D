reset
set si 0.7,0.7
set si rat 1.0
set term postscript eps size 5.5,5.5 enhanced color
set output "x-s.eps"
# set terminal pngcairo  background "#ffffff" fontscale 1.0 dashed size 640, 480 
# set output 'dashcolor.2.png'
set size square
set style function linespoints

set xtics font ", 24"
set ytics font ", 24"

set key font "Roman,30"
set xlabel "r/a_1" offset 0,-1 font "Roman,40"
set ylabel "{/Symbol s}" offset -1,-1 font "Roman,40"
set xrange[] noreverse nowriteback
set yrange[] noreverse nowriteback


#set ytics 0.0,0.2,1.2

set border 31 lw 2
set key spacing 6.5
#set key at 0.37,22.0,5.0

#set bmargin  6
#unset colorbox
i = 8
set grid

set lmargin at screen 0.11

set style line 1 lt 7 linecolor rgb "blue"
set style line 2 lt 2 lw 2 pt 3 ps 0.5 linecolor rgb "dark-green"

plot "numerical_x.dat" u 2:3 w l ls 1 lw 4 title '{/Symbol s_{11}}', "" u 2:7 w lp lt -1 pt 7 ps 1 title '{/Symbol s_{13}}'

#####################################################

set output "y-s.eps"

set xrange[] noreverse nowriteback
set yrange[-2.1:2.1] noreverse nowriteback

set border 31 lw 2
set key spacing 6.5
#set key at 0.40,9.5,5.0
set grid


plot "numerical_y.dat" u 2:3 w l ls 1 lw 4 title '{/Symbol s_{11}}', "" u 2:7 w lp lt -1 pt 7 ps 1 title '{/Symbol s_{13}}'

#####################################################

set output "z-s.eps"

set xrange[] noreverse nowriteback
set yrange[-2.2:1] noreverse nowriteback

set border 31 lw 2
set key spacing 6.5
#set key at 0.40,9.5,5.0
set grid


plot "numerical_z.dat" u 2:3 w l ls 1 lw 4 title '{/Symbol s_{11}}', "" u 2:7 w lp lt -1 pt 7 ps 1 title '{/Symbol s_{13}}'

#####################################################

set output "xy-s.eps"

set xrange[] noreverse nowriteback
set yrange[-2.2:1] noreverse nowriteback

set border 31 lw 2
set key spacing 6.5
set key at 3.9,-0.5,5.0
set grid


plot "numerical_xy.dat" u 2:3 w l ls 1 lw 4 title '{/Symbol s_{11}}', "" u 2:6 w lp lt 1 pt 6 ps 1 title '{/Symbol s_{23}}', "" u 2:7 w lp lt -1 pt 7 ps 1 title '{/Symbol s_{13}}'
