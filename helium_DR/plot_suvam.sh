




CMDDIR=~/calc/cmd

grep -v '0.000000e+00' <E_W_S_strontum.csv> outputdata2.csv

# Run FAC calculation
#python Kr_DR_He_plot.py  

# Turn csv file into plottable CS file
cat > gauein.dat << !
  E_W_S_strontum.csv
  Sigma.dat
       1                    IMODE
  0.0E-2  1.0E1   1.0E-4     EMIN, EMAX, DELE
  1.0E-4                    RTHRES
  5.0E-2                    HWIDTH
!
$CMDDIR/gauss.x
rm gauein.dat

# Plot CS with gnuplot
gnuplot  << !
set xrange[1:5]
set title "Total DR cross section"
set xlabel 'E [eV]'
set ylabel '{/Symbol s} [b]'
set key left top Left reverse
set key width -12
set term post enhanced "Times-Roman" 18
set out "Sigma-DR_suv.ps"
plot 'Sigma.dat' t "Sigma" w l lt 1 lw 3
!
ps2pdf Sigma-DR_suv.ps