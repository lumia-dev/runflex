#!/usr/bin/env bash

#SBATCH -t 12:00:00
#SBATCH -N 1
#SBATCH --output=$HOME/flexpart_%j.out

/proj/flexpart_aerosol/runflex/bin/runflex python $1 
