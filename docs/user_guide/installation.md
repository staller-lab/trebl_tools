# Installation

This analysis is designed to run on Savio to handle large data and computations. You can either run the code in a jupyter notebook on OOD or submit a job (see [Advanced Usage](advanced_usage.md) for submitting as a job).

## Option 1 (Savio): Use the shared lab environment

You can use the prebuilt lab conda environment directly:

```bash
conda activate /global/scratch/projects/fc_mvslab/conda/trebl_tools

# register a Jupyter kernel for this env (run once)
python -m ipykernel install --user --name trebl_tools_shared --display-name "trebl_tools (shared)"
```

## Option 2: Create your own conda environment

We still recommend creating your own conda environment for reproducibility and package isolation:

```bash
# clone the latest release 
git clone --branch v0.1.5 --depth 1 https://github.com/staller-lab/trebl_tools.git
cd trebl_tools

# create and activate conda env from the repo YAML
conda env create -f trebl_tools_env.yaml
conda activate trebl_tools_env  

# install the package (regular / non-editable)
pip install .

# install Jupyter kernel for this env
python -m ipykernel install --user --name trebl_tools_env --display-name "trebl_tools (v0.1.5)"
```

## Using Jupyter on Savio

1. Start a jupyter server session on Savio OOD

2. Open a notebook and select either:
   - `trebl_tools_shared` (shared lab environment), or
   - `trebl_tools_env` (your own environment)

3. The trebl_tools package will be available to import directly:

```python
import sys
import os
import glob

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import duckdb
from tqdm import tqdm

from trebl_tools import (
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
