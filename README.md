# TREBL Tools

Tools for TREBL analysis and barcode processing.

**For Staller Lab members on Savio.**

## Installation

### Recommended: Create Your Own Conda Environment

We recommend creating your own conda environment for better package management and easier Jupyter integration:

```bash
# 1. Clone the repository
git clone https://github.com/staller-lab/trebl_tools.git
cd trebl_tools

# 2. Create conda environment from the provided YAML file
conda env create -f trebl_tools_env.yaml
conda activate trebl_tools_env

# 3. Install trebl_tools package
pip install -e .

# 4. Install Jupyter kernel using your conda environment
python -m ipykernel install --user --name trebl_tools_env

# 5. Verify kernel is installed
jupyter kernelspec list