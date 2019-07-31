#!/bin/bash
#
#SBATCH --job-name=test
#SBATCH --output=res_18.txt
#
#SBATCH --ntasks=1
#SBATCH --time=1000:00
#SBATCH --mem-per-cpu=5000

./produceReco_18.sh
