#!/usr/bin/env bash
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH -n 32 

singularity run -H /home/x_guimo/PROJ/runflex/:/flexpart -B /proj/inversion/FLEXPART/meteo:/meteo -B /proj/inversion/FLEXPART/build/:/build -B /proj/inversion/FLEXPART/run:/run -B /proj/inversion/FLEXPART/results/:/output --app calcfootprints flexpart scripts/rc/singularity.rc flexpart_obslist.verify.csv
