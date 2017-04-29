#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=100gb,walltime=40:00:00
#PBS -M wrshoema@umail.iu.edu
#PBS -m abe
#PBS -j oe

#module load python
# keep it at 10 mutations in the population total per generation
#python /N/dc2/projects/Lennon_Sequences/2017_EvoDormReview/Python/runSimulations.py -6 -m -N 100 -M 100 -u 0.001 -G 100 -g 10000 -r 100
#python /N/dc2/projects/Lennon_Sequences/2017_EvoDormReview/Python/runSimulations.py -5 -m -N 100 -M 100 -u 0.001 -G 50 -g 1000 -r 100

python /Users/WRShoemaker/GitHub/EvoDormReview/Python/runSimulations.py -N 100 -M 100 -u 0.001 -G 50 -g 1000 -r 100
