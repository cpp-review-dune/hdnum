set term postscript eps enhanced color "Times-Roman,18"
#set style data linespoints
set key left top

set output "corona_daten.eps"
set grid
set xdata time
set timefmt "%d/%m/%y"
set xrange ["24/02/2020":"18/04/2020"]
set format x "%d.%m"
set xlabel "Datum"
set ylabel "Infizierte Personen"
plot "N.dat" u 1:2 lw 4 title "Infizierte"
! epstopdf corona_daten.eps

set output "corona_daten_logscale.eps"
set grid
set logscale y
set xdata time
set timefmt "%d/%m/%y"
set xrange ["24/02/2020":"18/04/2020"]
set format x "%d.%m"
set xlabel "Datum"
set ylabel "Infizierte Personen"
plot "N.dat" u 1:2 lw 4 title "Infizierte"
! epstopdf corona_daten_logscale.eps

set output "corona_vorhersage5_logscale.eps"
set grid
set logscale y
set xdata time
set timefmt "%d/%m/%y"
set xrange ["24/02/2020":"18/04/2020"]
set format x "%d.%m"
set xlabel "Datum"
set ylabel "Infizierte Personen"
plot "N.dat" u 1:2 lw 4 title "Daten", \
"fit1.dat" u 1:2 w l lw 4 title "Interpolation", \
"fit2.dat" u 1:2 w l lw 4 title "Ausgleichsgerade", \
"fit3.dat" u 1:2 w l lw 4 title "nichtl. Schaetzer"
! epstopdf corona_vorhersage5_logscale.eps

set output "corona_vorhersage5.eps"
set grid
set nologscale y
set xdata time
set timefmt "%d/%m/%y"
set xrange ["24/02/2020":"18/04/2020"]
set format x "%d.%m"
set xlabel "Datum"
set ylabel "Infizierte Personen"
plot "N.dat" u 1:2 lw 4 title "Daten", \
"fit1.dat" u 1:2 w l lw 4 title "Interpolation", \
"fit2.dat" u 1:2 w l lw 4 title "Ausgleichsgerade", \
"fit3.dat" u 1:2 w l lw 4 title "nichtl. Schaetzer"
! epstopdf corona_vorhersage5.eps

set output "corona_verdopplungszeiten.eps"
set grid
set nologscale y
set xdata time
set timefmt "%d/%m/%y"
set xrange ["24/02/2020":"18/04/2020"]
set format x "%d.%m"
set xlabel "Datum"
set ylabel "Verdopplungszeit [Tag]"
plot "verdopplungszeiten4.dat" u 1:2 w l lw 4 title "w=4", \
"verdopplungszeiten6.dat" u 1:2 w l lw 4 title "w=6",\
"verdopplungszeiten_naiv.dat" u 1:2 w l lw 4 title "naiv"
! epstopdf corona_verdopplungszeiten.eps
