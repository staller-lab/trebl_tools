TREBL Tools Documentation
=========================

Tools for TREBL (Transcriptional REporter Barcode Library) analysis and barcode processing.

**For Staller Lab members on Savio.**

Quick Links
-----------

* :doc:`user_guide/introduction` - Learn what TREBL-seq is
* :doc:`user_guide/installation` - Get started with installation
* :doc:`user_guide/analysis_setup` - Set up your analysis pipeline

Choose Your Analysis Path
-------------------------

Use the flowchart below to pick the right workflow before you start.

.. code-block:: text

   Do you have UMIs in your experiment?
   │
   ├── Yes ──► Is speed more important than maximum accuracy?
   │           │
   │           ├── Yes ──► Quick Start  (~2–4 hours, good for exploration/testing)
   │           │           • error_correction = False
   │           │           • umi_deduplication = "simple"
   │           │           Example: examples/notebooks/quick_start_example.ipynb
   │           │
   │           └── No  ──► Full Analysis  (~6–12 hours, publication-quality)
   │                       • error_correction = True
   │                       • umi_deduplication = "both"  (simple + directional)
   │                       Example: examples/notebooks/full_analysis_example.ipynb
   │
   └── No  ──► No-UMI Path  (see Advanced Usage)
               • AD_umi_object = None,  RT_umi_object = None
               • Set AD_reads_threshold / RT_reads_threshold manually

.. list-table:: At a glance: workflow comparison
   :header-rows: 1
   :widths: 30 23 23 24

   * - Parameter / Feature
     - Quick Start
     - Full Analysis
     - No-UMI
   * - ``error_correction``
     - ``False``
     - ``True``
     - ``False``
   * - ``umi_deduplication``
     - ``"simple"``
     - ``"both"``
     - *(omit UMI objects)*
   * - Typical runtime
     - ~2–4 hours
     - ~6–12 hours
     - varies
   * - Best for
     - Exploration / testing
     - Publication results
     - No UMI design
   * - Example notebook
     - `quick_start_example.ipynb <https://github.com/staller-lab/trebl_tools/blob/main/examples/notebooks/quick_start_example.ipynb>`_
     - `full_analysis_example.ipynb <https://github.com/staller-lab/trebl_tools/blob/main/examples/notebooks/full_analysis_example.ipynb>`_
     - :doc:`user_guide/advanced_usage`

User Guide
==========

.. toctree::
   :maxdepth: 2
   :caption: User Guide

   user_guide/introduction
   user_guide/installation
   user_guide/preprocessing
   user_guide/analysis_setup
   user_guide/step1
   user_guide/step2
   user_guide/trebl_experiment
   user_guide/advanced_usage

API Reference
=============

.. toctree::
   :maxdepth: 2
   :caption: API Documentation

   modules

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`