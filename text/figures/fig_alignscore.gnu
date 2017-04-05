load "include.gnu"
set terminal postscript portrait color dashed enhanced lw 1 size 16cm, 12cm font "Arial,9"

set output "fig_alignscore.eps"

set style data yerrorl

d = "../../data/"

set ytic 10000
set ylabel "Alignment score"
set xlabel "Sequence divergence"
c = 10
set xrange [0:.3]

set key samplen 6 spacing 1.3 vertical maxrows 3 at screen 0.9,0.99 width -3

# Outer multiplot margin. Units in fractions of drawing area.
ml=0.13                         # left
mr=0.96                         # right, counted from left, not closest edge
mb=0.1                         # bottom
mt=0.84                         # top, counted from bottom, not closest edge

# Spacing between panels
xspacing = 0.15
yspacing = 0.17

set multiplot \
    layout 2,2 rowsfirst \
    margins ml, mr, mb, mt \
    spacing xspacing, yspacing
# scale 1.2, 0.8   # not useful if setting margin and spacing


set label 2 "A" at graph -.2, 1.15 font "Arial-BoldMT,11"
set label 1 "Point mutations only" at graph 0.01,1.12
set ytic 5000
c = 16
plot[:][:] d."mutationrate_ri0.out" \
         u 1:c:c+1 ls 1 ti "Hirschberg",\
      "" u 1:c=c+2:c+1 ls 2 ti "Heuristic alignment",\
      "" u 1:c=c+2:c+1 ls 3 ti "Synapsing",\
      "" u 1:c=c+2:c+1 ls 4 ti "Synapsing (aligned)",\
      "" u 1:c=c+2:c+1 ls 5 ti "Heuristic 1-point",\
      "" u 1:c=c+2:c+1 ls 6 ti "Cut and splice",\
      "" u 1:c=c+2:c+1 ls 7 ti "SAGA"

unset key
set label 2 "B"
set label 1 "With indels"
c = 16
plot[:][:] d."mutationrate_ri0.05.out" \
         u 1:c:c+1 ls 1 ti "Hirschberg",\
      "" u 1:c=c+2:c+1 ls 2 ti "Heuristic alignment",\
      "" u 1:c=c+2:c+1 ls 3 ti "Synapsing",\
      "" u 1:c=c+2:c+1 ls 4 ti "Synapsing (aligned)",\
      "" u 1:c=c+2:c+1 ls 5 ti "Heuristic 1-point",\
      "" u 1:c=c+2:c+1 ls 6 ti "Cut and splice",\
      "" u 1:c=c+2:c+1 ls 7 ti "SAGA"

set ytic .1
set label 2 "C"
set label 1 "Point mutations only"
set ylabel "Crossover correctness"
c = 30
plot[:][:1.02] d."mutationrate_ri0.out" \
         u 1:c:c+1 ls 1 ti "Hirschberg",\
      "" u 1:c=c+2:c+1 ls 2 ti "Heuristic alignment",\
      "" u 1:c=c+2:c+1 ls 3 ti "Synapsing",\
      "" u 1:c=c+4:c+1 ls 5 ti "Heuristic 1-point",\
      "" u 1:c=c+4:c+1 ls 7 ti "SAGA"
#      "" u 1:c=c+2:c+1 ls 4 ti "Synapsing (aligned)",\
#      "" u 1:c=c+2:c+1 ls 6 ti "Cut and splice",\

set label 2 "D"
set label 1 "With indels"
c = 30
plot[:][:1.02] d."mutationrate_ri0.05.out" \
         u 1:c:c+1 ls 1 ti "Hirschberg",\
      "" u 1:c=c+2:c+1 ls 2 ti "Heuristic alignment",\
      "" u 1:c=c+2:c+1 ls 3 ti "Synapsing",\
      "" u 1:c=c+4:c+1 ls 5 ti "Heuristic 1-point",\
      "" u 1:c=c+4:c+1 ls 7 ti "SAGA"
#      "" u 1:c=c+2:c+1 ls 4 ti "Synapsing (aligned)",\
#      "" u 1:c=c+2:c+1 ls 6 ti "Cut and splice",\



unset multiplot

exit

set output "fig_alignscore_2.eps"
set multiplot layout 2,2



set ytic 5000
set ylabel "Alignment Score"
c = 8
set xlabel "Indel rate (mut = 0, {/Symbol a} = -1.5)"
plot[:][:] d."indelrate_m0.out" \
                        u 1:c:c+1 ti "Hirschberg",\
      "" u 1:c+2:c+3 ti "Heuristic",\
      "" u 1:c+4:c+5 ti "Synapsing"
set xlabel "Indel rate (mut = 0.05, {/Symbol a} = -1.5)"
plot[:][:] d."indelrate_m0.05.out" \
                        u 1:c:c+1 ti "Hirschberg",\
      "" u 1:c+2:c+3 ti "Heuristic",\
      "" u 1:c+4:c+5 ti "Synapsing"

set ylabel "Alignment Score 2"
c = 14
set xlabel "Indel rate (mut = 0, {/Symbol a} = -1.5)"
plot[:][:] d."indelrate_m0.out" \
                        u 1:c:c+1 ti "Hirschberg",\
      "" u 1:c+2:c+3 ti "Heuristic",\
      "" u 1:c+4:c+5 ti "Synapsing"
set xlabel "Indel rate (mut = 0.05, {/Symbol a} = -1.5)"
plot[:][:] d."indelrate_m0.05.out" \
                        u 1:c:c+1 ti "Hirschberg",\
      "" u 1:c+2:c+3 ti "Heuristic",\
      "" u 1:c+4:c+5 ti "Synapsing"
