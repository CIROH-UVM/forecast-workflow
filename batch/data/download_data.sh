#!/bin/bash

#SBATCH --partition=bluemoon
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=8G
#SBATCH --time=30:00:00
#SBATCH --job-name=<JOB_NAME>
#SBATCH --output=%x_%j.out
#SBATCH --mail-type=FAIL

set -x

source /users/n/b/nbeckage/miniconda3/etc/profile.d/conda.sh
conda activate forecast
export PYTHONPATH=~/ciroh/forecast-workflow

python download_<SOURCE>_<YEAR>.py
