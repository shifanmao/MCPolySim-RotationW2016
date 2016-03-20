#!/bin/bash
#SBATCH --job-name=KHOU
#SBATCH --output=wlcsim.out
#SBATCH --error=wlcsim.error
#SBATCH --time=48:00:00
#SBATCH --qos=normal
#SBATCH --mem=2000
#SBATCH --ntasks-per-node=16
./runwlcsim
