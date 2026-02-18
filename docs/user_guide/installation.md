# Installation

This analysis is designed to run on Savio to handle large data and computations. You can either run the code in a jupyter notebook on OOD or submit a job (see [Advanced Usage](advanced_usage.md) for submitting as a job).

## Recommended: Create Your Own Conda Environment

We recommend creating your own conda environment for better package management and easier Jupyter integration:

```bash
# 1. Clone the repository
git clone https://github.com/staller-lab/trebl_tools.git
cd trebl_tools

# 2. Create conda environment from the provided YAML file
# This will install all dependencies including Python, pandas, duckdb, umi_tools, etc.
conda env create -f trebl_tools_env.yaml
conda activate trebl_tools_env

# 3. Install trebl_tools package
pip install -e .

# 4. Install Jupyter kernel using your conda environment
python -m ipykernel install --user --name trebl_tools_env

# 5. Verify kernel is installed
jupyter kernelspec list
```

## Using Jupyter on Savio

1. Start a jupyter server session on Savio OOD

2. Open a notebook and select "trebl_tools_env" as the kernel

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
