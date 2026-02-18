# Installation

This analysis is designed to run on Savio to handle large data and computations. You can either run the code in a jupyter notebook on OOD or submit a job (see [Advanced Usage](advanced_usage.md) for submitting as a job).

## Recommended: Create Your Own Conda Environment

We recommend creating your own conda environment for better package management and easier Jupyter integration:

```bash
# 1. Create a new conda environment with Python 3.11
conda create -n trebl_env python=3.11
conda activate trebl_env

# 2. Install necessary conda packages (bioinformatics tools)
conda install -c bioconda umi_tools fastp

# 3. Install trebl_tools and its Python dependencies from git (v0.1.0-beta)
pip install git+https://github.com/staller-lab/trebl_tools.git@v0.1.0-beta

# 4. Install Jupyter kernel using your conda environment
python -m ipykernel install --user --name trebl_env

# 5. Verify kernel is installed
jupyter kernelspec list
```

## Using Jupyter on Savio

1. Start a jupyter server session on Savio OOD

2. Open a notebook and select "trebl_env" as the kernel

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
