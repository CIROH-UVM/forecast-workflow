#!/bin/bash

#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=24G
#SBATCH --time=18:00:00
#SBATCH --job-name=$job_name
#SBATCH --output=%x_%j.out
#SBATCH --mail-type=FAIL

set -xe

# Function to print error message on failure
failed_run() {
    echo "ERROR: RUN FAILED"
}

# Trap the ERR signal to call the error handler when a command fails
trap failed_run ERR

# source /users/p/c/pclemins/usr/local/miniforge3/etc/profile.d/conda.sh
source /users/n/b/nbeckage/miniconda3/etc/profile.d/conda.sh

# activate forecast environment
conda activate forecast

# export PYTHONPATH=/users/p/c/pclemins/repos/forecast-workflow
export PYTHONPATH=/users/n/b/nbeckage/ciroh/forecast-workflow

cd $run_dir

python -m models.aem3d.AEM3D_prep_worker --conf configuration.json
python -m models.aem3d.AEM3D_worker --conf configuration.json