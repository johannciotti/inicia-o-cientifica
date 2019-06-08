reset
set xrange [0.0:1.0]
set yrange [0:0.5]
set xtics  0.2
set ytics  0.2
set xlabel "{p}" offset 0,0.3 font 'Helvetica,35';
set ylabel "{g_{ex}} " offset 1.3 font 'Helvetica,35'
set cbrange [0.5:1]
set cbtics  0.1 offset -0.5,0.0 
set cblabel "{R} " offset -6.0,8.65 rotate by 0 left font 'Helvetica,35'
plot "POK_CV_F_N100_p_gexc_g1_gel0.1_-8.dat" u 1:2:4 with image t"
set style line 1 lc rgb "white" lt 1 lw 8 pt 7 ps 3.0
set terminal postscript eps enhanced color font 'Helvetica,25'
set output 'Fig6.eps'
replot "linha.dat" u 1:2  with linespoints ls 1 t"
