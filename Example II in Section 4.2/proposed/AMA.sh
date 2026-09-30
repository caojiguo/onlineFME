#!/bin/bash

Target_file1="./multi_logit.R"
LogDir="./logs"
mkdir -p "$LogDir" 
njob=8
Nrep=500
n_list=(10000 20000 30000 40000 50000)
ncores=1

export OMP_NUM_THREADS=$ncores
export OPENBLAS_NUM_THREADS=$ncores
export MKL_NUM_THREADS=$ncores
export NUMEXPR_NUM_THREADS=$ncores
export VECLIB_MAXIMUM_THREADS=$ncores

ccount=0

for n in "${n_list[@]}"
do
  for iloop in $(seq 1 $Nrep)
  do
      ccount=$((ccount+1))
      iflag=$((ccount % njob))

      if [ $iflag -eq 0 ]; then
          sleep 10s
          nohup R < "$Target_file1" --no-save --args $iloop $n \
          >> "$LogDir/Z1result_n${n}_rep_$(printf "%03d" $iloop).log" 2>&1
      else
          nohup R < "$Target_file1" --no-save --args $iloop $n \
          >> "$LogDir/Z1result_n${n}_rep_$(printf "%03d" $iloop).log" 2>&1 &
      fi
  done
done
