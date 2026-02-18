"""
TREBL Tools: Tools for TREBL analysis

Main classes for lab use:
- TreblPipeline: Run complete TREBL analysis workflows
"""

__version__ = "0.1.0"

# Main user-facing classes
from .pipelines import TreblPipeline

# Import submodules for user access
from . import preprocess
from . import finder
from . import initial_map
from . import map_refiner
from . import complexity
from . import error_correct
from . import plotting
from . import umi_deduplicate

__all__ = [
    "TreblPipeline",
    "preprocess",
    "finder",
    "initial_map",
    "map_refiner",
    "complexity",
    "error_correct",
    "plotting",
    "umi_deduplicate",
]
