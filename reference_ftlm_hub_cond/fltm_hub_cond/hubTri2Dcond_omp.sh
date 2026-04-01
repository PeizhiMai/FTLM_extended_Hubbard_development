#!/bin/sh                                                                                                                             

# compiling

\rm hubTri2Dcond_omp.e
echo compiling...

rm report
echo "host:" >> report
hostname >> report

ifort -O3 -qopenmp -mcmodel=large -no-wrap-margin -heap-arrays  hubTri2Dcond_omp.f  -o hubTri2Dcond_omp.e >> report

# running

hostname
date

nthreads=2
export OMP_NUM_THREADS=$nthreads
echo Number of threads $OMP_NUM_THREADS

echo running...

ulimit -s unlimited

time ./hubTri2Dcond_omp.e >& hubTri2Dcond_omp.log
     
\rm hubTri2Dcond_omp.e

