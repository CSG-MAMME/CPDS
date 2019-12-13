#!/bin/bash
OMP_NUM_THREADS=$1 ./heatomp test.dat #| tee jacobi_out_$1.txt
