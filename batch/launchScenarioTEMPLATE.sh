#!/bin/bash

#SBATCH --partition=bluemoon
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=8G
#SBATCH --time=30:00:00
#SBATCH --job-name=<YOUR_JOB_NAME>
#SBATCH --output=%x_%j.out
#SBATCH --mail-type=ALL

source /users/n/b/nbeckage/miniconda3/etc/profile.d/conda.sh
# source /users/p/c/pclemins/usr/local/miniforge3/etc/profile.d/conda.sh
conda activate forecast

# change to whatever directory you want to put your scenario runs in
cd /netfiles/ciroh/7dayHABsHindcast/<FEE_VERSION>/<CQ_PARADIGM>/<SCENARIO_DIR>/

# note that you must have an experiment-specific configuration file in your scenario directory
python ../../launchExperiment.py experiment_config.json
