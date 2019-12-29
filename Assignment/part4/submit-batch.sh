#!/bin/bash
# @ job_name		= heatCUDA.extrae
# @ partition		= debug
# @ initialdir		= .
# @ output		= heatCUDA.extrae.%j.out
# @ error		= heatCUDA.extrae.%j.err
# @ total_tasks		= 1
# @ gpus_per_node	= 1
# @ wall_clock_limit	= 00:02:00

executable=./heatCUDA
# Define threads per block
txb=1
max_txb=16

# module load extrae

while (test $txb -le $max_txb)
  do
    echo "Running with $txb threads per block."
    ./trace.sh ${executable} test.dat -t $txb
    txb=`expr $txb \* 2`
  done
