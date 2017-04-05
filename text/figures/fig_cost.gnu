load "include.gnu"
set terminal postscript portrait color enhanced lw 1 size 16cm, 12cm font "Arial,9"

set output "fig_cost.eps"

# Outer multiplot margin. Units in fractions of drawing area.
ml=0.1                         # left
mr=0.96                         # right, counted from left, not closest edge
mb=0.12                         # bottom
mt=0.85                         # top, counted from bottom, not closest edge

# Spacing between panels
xspacing = 0.18
yspacing = 0.17

set multiplot \
    layout 2,2 rowsfirst \
    margins ml, mr, mb, mt \
    spacing xspacing, yspacing
# scale 1.2, 0.8   # not useful if setting margin and spacing

set ylabel "-log(-fitness)"
set xlabel "Generations"

# rätta till linjetyper. Funkar inte bra med for-loopen

set key samplen 6 spacing 1.3 vertical maxrows 4 at graph 1,1.47 width -3

dir = "../../evodata2/"
d = 3*3-1;
cs = "0 0.1 0.3 0.6"
cs = "0 0.3"
set ytic auto
set samp 10

set xtic 1e4
set label 2 "A" at graph -.2, 1.25 font "Arial-BoldMT,11"
plot [][-3.3:-.9] \
     dir."0/evolution_average_m1.out"   u 1:(-log10(-column(d))) w l dt 2 lw 2 lc rgb "black" ti "No crossover", \
     dir."0.3/evolution_average_m1.out" u 1:(-log10(-column(d))) w l ls 2 ti "Heuristic alignment", \
     dir."0.3/evolution_average_m2.out" u 1:(-log10(-column(d))) w p pt 3 lc rgbcolor synaps ps 0.3 ti "Synapsing"
#     -log10(30) w l dt 3 lw 2 lc rgb "orange" ti "Target, fitness = -30"

#set ylabel "Fitness"
#cs = "0 0.05 0.1 0.3 0.6 0.7"
#set ytics auto
#plot [15e3:20e3][:] for [c in cs] dir.c."/evolution_average_m1.out" u 1:((column(d))) w l lt c*5 lc rgb "black" ti c, \
#	for [c in cs] dir.c."/evolution_average_m2.out" u 1:((column(d))) w p lt c*5 lc rgb "black" pt 5+c*10 ps .3 noti

set key at graph 1,1.39


f=-30
set label 2 "B"
set ylabel "Generations to target"
set xlabel "Crossover probability"
set ytic 6000
num="x{0,0.04,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8}"
set xtic auto
plot [-.01:.81] \
        "<cat ../../evodata/generations_f".f."_m1_".num.".out |sort -g" u 1:3:4 ls 2 w yerrorl ti "Heuristic alignment", \
        "<cat ../../evodata/generations_f".f."_m2_".num.".out |sort -g" u 1:3:4 ls 3 w yerrorl ti "Synapsing"

unset key
set label 2 "C"
set ytic 100
set ylabel "CPU time to target (s)"
plot [-.01:.81] \
        "<cat ../../evodata/generations_f".f."_m1_".num.".out |sort -g" u 1:18:19 ls 2 w yerrorl ti "Heuristic alignment", \
        "<cat ../../evodata/generations_f".f."_m2_".num.".out |sort -g" u 1:18:19 ls 3 w yerrorl ti "Synapsing

set label 2 "D"
set ytic 50
set ylabel "CPU time, hypothetical (s)"
s=0.015
plot [-.01:.51] \
        "<cat ../../evodata/generations_f".f."_m1_".num.".out |sort -g" u 1:($18+s*$3):($19+s*$4) ls 2 w yerrorl ti "Heuristic alignment", \
        "<cat ../../evodata/generations_f".f."_m2_".num.".out |sort -g" u 1:($18+s*$3):($19+s*$4) ls 3 w yerrorl ti "Synapsing
