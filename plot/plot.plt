set terminal aqua font "Helvetica,16" enhanced solid'
set logs
set yr[1.e-32:1.e0]
set xr[0.5e-3:0.5]
#set format y "%1.0e"
set format y "10^{%L}"
set key bottom right font "Helvetica,12" spacing 1.0
set grid
set title "8-figure, 10 orbits (T = 2{/Symbol p})"
set xlabel "{/Symbol D}t"
set ylabel "max |{/Symbol D}E/E|"

pl \
  '4th-pec1.dat' w linesp  lt 7 lw 1 pt  8 ti "4th PEC", \
  '4th-pec2.dat' w linesp  lt 7 lw 2 pt 11 ti "4th P(EC)^2", \
  '6th-pec1.dat' w linesp  lt 2 lw 1 pt  8 ti "6th PEC", \
  '6th-pec2.dat' w linesp  lt 2 lw 2 pt 11 ti "6th P(EC)^2", \
  '8th-pec1.dat' w linesp  lt 3 lw 1 pt  8 ti "8th PEC", \
  '8th-pec2.dat' w linesp  lt 3 lw 2 pt 11 ti "8th P(EC)^2", \
  '10th-pec1.dat' w linesp lt 8 lw 1 pt  8 ti "10th PEC", \
  '10th-pec2.dat' w linesp lt 8 lw 2 pt 11 ti "10th P(EC)^2", \
  '12th-pec1.dat' w linesp lt 4 lw 1 pt  8 ti "12th PEC", \
  '12th-pec2.dat' w linesp lt 4 lw 2 pt 11 ti "12th P(EC)^2", \
  '14th-pec1.dat' w linesp lt 1 lw 1 pt  8 ti "14th PEC", \
  '14th-pec2.dat' w linesp lt 1 lw 2 pt 11 ti "14th P(EC)^2" , \
  '16th-pec1.dat' w linesp lt 0 lw 1 pt  8 ti "16th PEC", \
  '16th-pec2.dat' w linesp lt 0 lw 2 pt 11 ti "16th P(EC)^2", \
  '16th-pec3.dat' w linesp lt 0 lw 2 pt  3 ti "16th P(EC)^3" #, \
