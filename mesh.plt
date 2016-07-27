set terminal wxt
set xyplane 0
set view equal xyz
unset border
unset xtics
unset ytics
unset ztics
set lmargin 0
set rmargin 0
set tmargin 0
set bmargin 0

splot 'primalgrid.dat' using 1:2:3:($4-$1):($5-$2):($6-$3) with vectors nohead notitle
