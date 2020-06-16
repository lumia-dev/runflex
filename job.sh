#!/usr/bin/env bash

#SBATCH -t 12:00:00
#SBATCH -N 1
#SBATCH -n 25
#SBATCH -A snic2019-7-72
#SBATCH --output=$HOME/flexpart_%j.out

/proj/flexpart_aerosol/runflex_moa/runflex/bin/runflex python  /proj/flexpart_aerosol/runflex_moa/runflex/scripts/calcFootprints.py 
