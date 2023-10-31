# forecast-workflow
CIROH @ UVM Forecast Workflow Model and Data Preparation Components

## Repository workflow for developer's reference
### 1. `models/aem3d/AEM3D_prep_worker.py`
This is the file that loads `models/aem3d/AEM3D_prep_IAM.py` as a module, which sets up logging and then runs `AEM3D_prep_IAM()` to setup data and get all inputs ready for the model run. Note that `AEM3D_prep_IAM()` *does not* download any data by default.
### 2. `models/aem3d/AEM3D_worker.py`
This file runs the AEM3D model with logging. See `main()` to see exactly what is going on, but so far it looks like the only settings being directly passed to `aem3d_openmp` besides an empty list `args` is `_iter=True`. I'm not sure how the model is inheriting all of the setup done in step 1. 