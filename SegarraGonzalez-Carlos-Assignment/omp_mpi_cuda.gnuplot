#!/usr/bin/gnuplot
#load "styles.inc"
set terminal postscript color eps enhanced font 22 #size 5,2
set output 'omp_mpi_cuda.eps'
set datafile separator ","
set bmargin 4
set tmargin 2.5
set xrange [-10:200]
set xtics 20
set xlabel "Delay to Timeout Ratio [\%]"
#set xtics rotate by -45 offset 0.0, 0.0 font ",12"
set yrange [0:3.5]
set y2range [0:3.5]
set ytics 1
set ylabel "Average Execution Time [s]" offset -1.0, -1.5 
set ytics nomirror
set grid y
set key left

# Line Styles: check line and point types in PostScript terminal here
# https://stackoverflow.com/questions/19412382/gnuplot-line-types
set style line 01 lt 1 lw 3 lc rgb "#fd9103" pt 4
set style line 11 lt 1 lw 3 dt 2 lc rgb "#FD9103" pt 4
set style line 02 lt 1 lw 3 lc rgb "#8803FD" pt 8
set style line 12 lt 2 lw 3 dt 2 lc rgb "#8803FD" pt 8
set style line 03 lt 1 lw 3 lc rgb "#BC1B36" pt 19
set style line 13 lt 2 lw 3 dt 2 lc rgb "#BC1B36" pt 19

#set title "{/bold Effect of Including Delays in the Acceptor}" offset 0.0,-0.5
set multiplot layout 2,1 \
              margins 0.1, 0.9, 0.15, 0.9 \
              spacing 0, 0.15
    
# OMP
set xrange [1:16]
set xtics (1, 2, 4, 8, 16)
set logscale x
set xlabel "# OMP\\_THREADS" font ",18"
set yrange [0:6]
set ytics 2
set ylabel "Execution Time [s]" offset 1,-6
set y2range [0:8]
set y2tics 2
set key at screen 0.0, screen 0.99 vertical maxrows 1 width -2 #font ",18"
set y2label "Speed Up" offset -1.8,-5.0 rotate by -90
plot './omp_data.dat' using 1:2 w linespoints ls 01 title 'Jacobi' axis x1y1, \
     './omp_data.dat' using 1:3 w linespoints ls 11 notitle '' axis x1y2, \
     './omp_data.dat' using 1:4 w linespoints ls 02 title 'Red Black' axis x1y1, \
     './omp_data.dat' using 1:5 w linespoints ls 12 notitle '' axis x1y2, \
     './omp_data.dat' using 1:6 w linespoints ls 03 title 'Gauss Seidel' axis x1y1, \
     './omp_data.dat' using 1:7 w linespoints ls 13 notitle '' axis x1y2

unset key
set xrange [1:16]
set logscale x
set ylabel ""
set y2label ""
set xlabel "# MPI Processes"
set yrange [0:6]
set y2range [0:8]
plot './mpi_data.dat' using 1:2 w linespoints ls 01 title 'Jacobi' axis x1y1, \
     './mpi_data.dat' using 1:3 w linespoints ls 11 notitle '' axis x1y2, \
     './mpi_data.dat' using 1:4 w linespoints ls 03 title 'Gauss' axis x1y1, \
     './mpi_data.dat' using 1:5 w linespoints ls 13 notitle '' axis x1y2

!epstopdf 'omp_mpi_cuda.eps'
!rm 'omp_mpi_cuda.eps'
