"""
TREBL Tools: Tools for TREBL analysis 

Main classes for lab use:
- TreblPipeline: Run complete TREBL analysis workflows
"""

__version__ = "0.1.0"

# Main user-facing classes
from .pipelines import TreblPipeline

__all__ = [
    "TreblPipeline",
]