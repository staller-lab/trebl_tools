# TREBL Tools

Tools for TREBL analysis and barcode processing.

**For Staller Lab members on Savio.**

## Installation

### Recommended: Create Your Own Conda Environment

We recommend creating your own conda environment for better package management and easier Jupyter integration:

```bash
# 1. Create a new conda environment with Python 3.11
conda create -n trebl_env python=3.11
conda activate trebl_env

# 2. Clone the repository
git clone https://github.com/staller-lab/trebl_tools.git
cd trebl_tools

# 3. Install trebl_tools and its Python dependencies
pip install -e .

# 4. Install Jupyter kernel using your conda environment
python -m ipykernel install --user --name trebl_env

# 5. Verify it's installed
jupyter kernelspec list