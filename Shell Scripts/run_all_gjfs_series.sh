#!/bin/bash
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=24:00:00
#SBATCH --job-name=m062x
#SBATCH --account=ucb273_peak1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=userid@colorado.edu

module purge
module load "gaussian/16_avx2"

FILES=*.gjf
for f in $FILES
do
 name=$(echo "$f" | cut -f 1 -d '.')
 g16 <$name.gjf> $name.log
  echo "Processing $name file..."
done

date
