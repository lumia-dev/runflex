#!/usr/bin/env bash

#SBATCH -t 02:00:00
#SBATCH -J flexpart_test
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -A snic2019-7-72
#SBATCH -o=$HOME/flexpart_%j.out
#SBATCH -e=$HOME/flexpart_%j.err

/proj/flexpart_aerosol/runflex_moa/runflex/bin/runflex python  /proj/flexpart_aerosol/runflex_moa/runflex/scripts/calcFootprints.py
