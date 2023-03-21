#!/bin/bash
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=24:00:00
#SBATCH --job-name=Rotors
#SBATCH --account=ucb273_peak1
#SBATCH --array=1-36
#SBATCH --mail-type=ALL
#SBATCH --mail-user=userid@colorado.edu


module purge
module load "gaussian/16_avx2"

FILES=$(ls *.gjf | sed -n ${SLURM_ARRAY_TASK_ID}p)
echo "Processing $FILES ...."

g16 $FILES

echo "This run performed on:"
date