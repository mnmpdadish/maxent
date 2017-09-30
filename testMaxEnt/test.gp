
set terminal wxt size 1400,800

set multiplot

unset logscale
set ylabel 'A(w)'
set ylabel 'w'
set tmargin at screen 0.1; set lmargin at screen 0.1; set rmargin at screen 0.95; set bmargin at screen 0.95;
plot '-' u 1:2 w lp t ''
1 2
2 2.2
3 2.1
4 1.7
e

unset yrange
set format y ''; unset ylabel
set format x ''; unset xlabel
set tmargin at screen 0.6; set lmargin at screen 0.7; set rmargin at screen 0.95; set bmargin at screen 0.95;
set logscale xy
plot '-' u 1:2 w lp t '' ps 0
1 2
2 2.2
3 2.1
4 1.7
e



unset multiplot

pause -1
