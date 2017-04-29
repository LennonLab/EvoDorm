#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=100gb,walltime=30:00:00
#PBS -M wrshoema@umail.iu.edu
#PBS -m abe
#PBS -j oe

module load python

python /Users/WRShoemaker/GitHub/EvoDormReview/Python/runSimulations.py -6 -s -N 100 -M 10 -u 0.0001 -G 100 -g 100 -c 0 -r 1000
