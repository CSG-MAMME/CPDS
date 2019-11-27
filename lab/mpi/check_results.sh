#!/bin/bash
solvers="jacobi red_black"

for sol in $solvers;
do
    echo "---------- Checking $sol solver! ----------"
    echo "Difference between output images, running diff:"
    diff ./initial-results/heat_$sol.ppm ./omp-par/heat_$sol.ppm
done
