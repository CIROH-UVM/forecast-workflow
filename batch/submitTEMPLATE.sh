#!/bin/bash

#SBATCH --partition=bluemoon
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=8G
#SBATCH --time=18:00:00
#SBATCH --job-name=$job_name
#SBATCH --output=%x_%j.out
#SBATCH --mail-type=FAIL

set -x

source /users/n/b/nbeckage/miniconda3/etc/profile.d/conda.sh
conda activate forecast
# export PYTHONPATH=/users/p/c/pclemins/repos/forecast-workflow
export PYTHONPATH=/users/n/b/nbeckage/ciroh/forecast-workflow

cd $run_dir

python -m models.aem3d.AEM3D_prep_worker --conf configuration.json
python -m models.aem3d.AEM3D_worker --conf configuration.json