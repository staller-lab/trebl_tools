# Installation

This analysis is designed to run on Savio to handle large data and computations. You can either run the code in a jupyter notebook on OOD or submit a job (see [Advanced Usage](advanced_usage.md) for submitting as a job).

## Setting up Jupyter on Savio

1. Start a jupyter server session on Savio.

2. From the terminal, run this command:

```bash
python -m ipykernel install --user --name trebl_env --prefix /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/conda/trebl_env
```

3. Then, open a notebook and select "trebl_env" as the kernel.

## Required Imports

Copy and paste the following imports in the first cell:

```python
import sys
import os
import glob

# Points to conda python path for correct package installs
sys.path.insert(0, "/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/conda/trebl_env/lib/python3.11/site-packages")

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import duckdb
from tqdm import tqdm

# Points to location of scripts
sys.path.append("/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/TREBL/")
from scripts import (
    initial_map,
    map_refiner,
    complexity,
    finder,
    preprocess,
    error_correct,
    plotting,
    umi_deduplicate,
    pipelines
)
```

**Note:** The `.yaml` file corresponding to the conda environment can be found in the repository.
