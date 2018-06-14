set encoding iso_8859_1
set terminal postscript enhanced color  "Helvetica" 16
set output "plot05.ps"

set size square 1,1
set xzeroaxis

set style line  1  lt  3         lc rgb "blue"     lw  4
set style line  2  lt  6         lc rgb "red"      lw  8
set style line  3  lt  9         lc rgb "green"    lw 12
set style line  4  lt  4         lc rgb "green"    lw  4
set style line  5  lt  5         lc rgb "magenta"  lw  4
set style line  6  lt  6         lc rgb "orange"   lw  4
set style line  7  lt  7         lc rgb "cyan"     lw  4
set style line  8  lt  8         lc rgb "coral"    lw  4
set style line  9  lt  9         lc rgb "brown"    lw  4
set style line 10  lt 10  pt  1  lc rgb "black"    lw  4
set style line 11  lt 11  pt  2  lc rgb "blue"     lw  4
set style line 12  lt 12  pt  3  lc rgb "red"      lw  4
set style line 13  lt 13  pt  4  lc rgb "green"    lw  4
set style line 14  lt 14  pt  5  lc rgb "magenta"  lw  4
set style line 15  lt 15  pt  6  lc rgb "orange"   lw  4
set style line 16  lt 16  pt  7  lc rgb "cyan"     lw  4
set style line 17  lt 17  pt  8  lc rgb "coral"    lw  4
set style line 18  lt 18  pt  9  lc rgb "brown"    lw  4
set style line 19  lt 19  pt 10  lc rgb "black"    lw  4
set style line 20  lt 20  pt 11  lc rgb "blue"     lw  4

set tics scale 1.0
set ytics border nomirror norotate 10000

set mxtics 2
set mytics 2

set xrange [180:250]
set yrange [-18000:24000]

set title "CD Spectra of Proteins"
show label
set xlabel "Wavelength" 0,-1 font "Helvetica, 22"
set ylabel "Intensity" 0,-1 font "Helvetica, 22"

set key right top outside Left reverse samplen 4 spacing 5 title "" nobox

plot "1rhd.exp.cd"    title "helical proteins"    with lines ls  1 , \
     "4gcr.exp.cd"    title "beta sheet proteins" with lines ls  2 , \
     "2cga.exp.cd"    title "PPII proteins"       with lines ls  3 