#!/bin/sh

nthreads=4

prog=cond_spect_omp

fname=$prog.f      # fortran file name
exename=$prog.e    # Executable name
repname=report1    # Reports name  
logname=$prog.log  # Log name 

# compiling 
echo compiling to $exename...
\rm $exename

#module load intel/2017
ifort -O3 -qopenmp -mcmodel=large -no-wrap-margin -heap-arrays $fname -o $exename >& $repname

hostname 
date

echo File name: $fname 
echo Log. file name: $logname 
echo Reserved cpus: $nthreads 

export OMP_NUM_THREADS=$nthreads 

ulimit -s unlimited 

#module load intel/2017
time ./$exename >& $logname

\rm $exename 

