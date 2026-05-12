# TREBL Tools

Tools for TREBL analysis and barcode processing.

**For Staller Lab members on Savio.**

## Installation

### Option 1 (Savio): Use the shared lab environment

You can use the prebuilt lab conda environment directly:

```bash
conda activate /global/scratch/projects/fc_mvslab/conda/trebl_tools

# register a Jupyter kernel for this env (run once)
python -m ipykernel install --user --name trebl_tools_shared --display-name "trebl_tools (shared)"
```

### Option 2: Create your own conda environment

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
