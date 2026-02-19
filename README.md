# TREBL Tools

Tools for TREBL analysis and barcode processing.

**For Staller Lab members on Savio.**

## Installation

### Recommended: Create Your Own Conda Environment

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