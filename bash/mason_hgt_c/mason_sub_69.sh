#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=100gb,walltime=30:00:00
#PBS -M wrshoema@umail.iu.edu
#PBS -m abe
#PBS -j oe

module load python

python /Users/WRShoemaker/GitHub/EvoDormReview/Python/runSimulations.py -5 -s -N 100 -M 100 -u 0.001 -G 50 -g 10 -c 69 -r 1
