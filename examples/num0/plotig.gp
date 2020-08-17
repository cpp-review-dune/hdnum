set term postscript eps enhanced color "Times-Roman,18"
#set style data linespoints
set key left top

set output "ig_rechteseite.eps"
set grid
set xlabel "x"
set ylabel "f(x)"
plot "rechteseite.dat" u 2:3 w l lw 4 title "f(x)"
! epstopdf ig_rechteseite.eps

set output "ig_loesungen.eps"
set grid
set xlabel "x"
set ylabel "u(x)"
plot "loesung_350.dat" u 2:3 w l lw 4 title "gamma=0.35",\
"loesung_500.dat" u 2:3 w l lw 4 title "gamma=0.5",\
"loesung_750.dat" u 2:3 w l lw 4 title "gamma=0.75"
! epstopdf ig_loesungen.eps

