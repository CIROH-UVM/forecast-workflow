#!/bin/bash

# Define the parent directory - the folder in your current working directory you want to comb through
# usually a year's worth of model runs: 2023/, 2022, etc
parent_dir="year"
# Define job prefix - defined in the scenario's experiment_config.json
job_pref="QM#_"
# Define error message - last line of your model run's out file that identifies failed runs
# note that carrot only looks for string at the beginning of the line (regex)
error="Exception: AEM3D failure."

# Loop through each subdirectory
for dir in "$parent_dir"/*; do
    if [ -d "$dir" ]; then
        # Check if the subdirectory contains the specified file
        if [ -e "$dir"/"$job_pref"* ]; then
            # Get the last line of the file and check if they contain the specified text
            if tail -n 1 "$dir"/"$job_pref"* | grep -q "$error"; then
                # If the text is found, remove the entire directory
                echo "Removing directory: $dir"
				tail -n 1 "$dir"/"$job_pref"*
				# CAUTION - uncomment to actually delete entire run dir - THIS CANNOT BE UNDONE!
                # rm -rf "$dir"
            fi
        fi
    fi
done
