#!/bin/bash
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=24:00:00
#SBATCH --job-name=MESS 
#SBATCH --account=ucb273_peak1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=userid@colorado.edu

module purge
module use /projects/nila3952/public/modules
module load miniconda
source activate mess-env
mess M1M_MESS_File_v8-0-6.inp

echo "This run performed on:"
date
