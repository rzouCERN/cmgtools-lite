#!/bin/bash
#
#SBATCH --job-name=test
#SBATCH --output=res_17.txt
#
#SBATCH --ntasks=1
#SBATCH --time=100:00:00
#SBATCH --mem-per-cpu=1000

./produceReco_17.sh
