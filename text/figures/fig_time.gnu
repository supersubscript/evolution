load "include.gnu"
set terminal postscript portrait color enhanced lw 1 size 16cm, 14cm font "Arial,9"

set output "fig_time.eps"

set style data yerrorl
set key at graph 0.65, .8

d = "../../data/"

#set xtic
unset mxtics

set key samplen 7 spacing 1.4 vertical maxrows 3 at screen 0.8,0.99 reverse Left

# Outer multiplot margin. Units in fractions of drawing area.
ml=0.1                         # left
mr=0.96                         # right, counted from left, not closest edge
mb=0.1                         # bottom
mt=0.83                         # top, counted from bottom, not closest edge

# Spacing between panels
xspacing = 0.15
yspacing = 0.15

set multiplot \
    layout 2,2 rowsfirst \
    margins ml, mr, mb, mt \
    spacing xspacing, yspacing
# scale 1.2, 0.8   # not useful if setting margin and spacing

set samp 5

set format y "10^{%T}"
set format x "10^{%T}"
set logs
set label 2 "A" at graph -.2, 1.15 font "Arial-BoldMT,11"
set label 1 "Point mutations only" at graph 0.01,1.12
set xlabel "Genome length"
set ylabel "CPU time (s)" offset 0.,0
c = 2
plot[1e3:1e6][:] d."lengths_m0.1_i0.out" u 1:c:c+1 ti "Hirschberg " ls 1,\
      d."lengths_m0.05_i0.out" u 1:c+2:c+3 ti "Heuristic alignment" ls 2,\
      "" u 1:c+4:c+5 ti "Synapsing" ls 3,\
      9e-9*x*x ti "O(N^2)" ls 10 ,\
      4.5e-7*x*log(x) ti "O(N ln N)" ls 11

unset key

set label 2 "B"
set label 1 "With insertions/deletions"
set xlabel "Genome length"
plot[1e3:1e6][:] d."lengths_m0.1_i0.005.out" u 1:c:c+1 ti "Hirschberg" ls 1,\
      d."lengths_m0.05_i0.0025.out" u 1:c+2:c+3 ti "Heuristic alignment" ls 2, \
      "" u 1:c+4:c+5 ti "Synapsing" ls 3,\
      .7e-9*x*x*log(x) ti "O(N^2)" ls 10,\
      4.5e-7*x*log(x) ti "O(N ln N)" ls 11

#unset format
unset logs

set logs y


set label 2 "C"
#set format x "%.0f%%"
unset format x
set ylabel "CPU time (s)" offset 0,0
set ytic auto
set label 1 "Point mutations only"
set xtic auto
set xlabel "Sequence divergence"
c = 2
plot[:][3e-3:3e-1] d."mutationrate_ri0.out" \
         u 1:c:c+1 ls 1 ti "Hirschberg",\
      "" u 1:c=c+2:c+1 ls 2 ti "Heuristic alignment",\
      "" u 1:c=c+2:c+1 ls 3 ti "Synapsing",\
      "" u 1:c=c+2:c+1 ls 4 ti "Synapsing (aligned)",\
      "" u 1:c=c+2:c+1 ls 5 ti "Heuristic 1-point",\
      "" u 1:c=c+2:c+1 ls 6 ti "Cut and splice",\
      "" u 1:c=c+2:c+1 ls 7 ti "SAGA"

set label 2 "D"
set label 1 "With insertions/deletions"
set xlabel "Sequence divergence"
c = 2
plot[:][3e-3:3e-1] d."mutationrate_ri0.05.out" \
         u 1:c:c+1 ls 1 ti "Hirschberg",\
      "" u 1:c=c+2:c+1 ls 2 ti "Heuristic alignment",\
      "" u 1:c=c+2:c+1 ls 3 ti "Synapsing",\
      "" u 1:c=c+2:c+1 ls 4 ti "Synapsing (aligned)",\
      "" u 1:c=c+2:c+1 ls 5 ti "Heuristic 1-point",\
      "" u 1:c=c+2:c+1 ls 6 ti "Cut and splice",\
      "" u 1:c=c+2:c+1 ls 7 ti "SAGA"
