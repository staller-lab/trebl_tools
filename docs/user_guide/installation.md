# Installation

This analysis is designed to run on Savio to handle large data and computations. You can either run the code in a jupyter notebook on OOD or submit a job (see [Advanced Usage](advanced_usage.md) for submitting as a job).

## Recommended: Create Your Own Conda Environment

We recommend creating your own conda environment for better package management and easier Jupyter integration:

```bash
# clone the latest release 
git clone --branch v0.1.0 --depth 1 https://github.com/staller-lab/trebl_tools.git
cd trebl_tools

# create and activate conda env from the repo YAML
conda env create -f trebl_tools_env.yaml
conda activate trebl_tools_env  

# install the package (regular / non-editable)
pip install .

# install Jupyter kernel for this env
python -m ipykernel install --user --name trebl_tools_env --display-name "trebl_tools (v0.1.0)"
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
