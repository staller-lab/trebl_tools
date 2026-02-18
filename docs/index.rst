TREBL Tools Documentation
=========================

Tools for TREBL analysis.

**For Staller Lab members on Savio.**

Installation
------------

.. code-block:: bash

   git clone https://github.com/staller-lab/trebl_tools.git
   cd trebl_tools
   pip install --user -e .

Quick Start
-----------

.. code-block:: python

   from trebl_tools import Barcode, TreblPipeline

   # Define barcodes
   ad_bc = Barcode(name="AD_BC", preceder="ATCG", post="GCTA", length=20)

   # Run pipeline
   pipeline = TreblPipeline(
       db_path="experiment.duckdb",
       output_dir="results/"
   )

   pipeline.trebl_experiment_analysis(
       AD_seq_files=["AD.fastq"],
       REP_seq_files=["REP.fastq"],
       AD_bc_objects=[ad_bc],
       REP_bc_objects=[ad_bc],
       design_file_path="design.csv"
   )

API Reference
=============

.. toctree::
   :maxdepth: 2

   modules

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`