#!/bin/bash
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=24:00:00
#SBATCH --job-name=kinbot_isooctane_P9_oxi
#SBATCH --account=ucb273_peak1

dos2unix kinbot.json
dos2unix job_template.text
dos2unix Transfer_Files.sh
chmod 777 Transfer_Files.sh

module purge
module use /projects/nila3952/public/modules
module load miniconda
module load gaussian/16_avx2
source activate Kinbot_Env

kinbot kinbot.json 

./Transfer_Files.sh