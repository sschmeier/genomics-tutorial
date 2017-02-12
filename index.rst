.. Genomics Tutorial documentation master file, created by
   sphinx-quickstart on Sat Nov 19 11:54:02 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Genomics Tutorial
=================

This is an introductory tutorial for learning genomics mostly on the command-line.
You will learn how to analyse next-generation sequencing (NGS) data.
The data you will be using is actual research data.
The final aim is to identify the genome variations in evolved lines of wild yeast that can explain the observed biological phenotypes.

The tutorial workflow is summarised in :numref:`fig-workflow`.

.. _fig-workflow:
.. figure:: images/workflow.png

   The tutorial will follow this workflow.



During this tutorial you will learn to:

- Make use of a UNIX-based computer environment
- Check the data quality of an NGS experiment
- Create a genome assembly of the ancestor based on NGS data
- Annotate a newly derived genome
- Map NGS reads of evolved lines to the ancestral reference genome
- Call genome variations/mutations in the evolved lines
- Identify the genes responsible for the observed evolved phenotypes

.. .. only:: builder_html

..   A printable PDF version of this tutorial can be downloaded :download:`here <_static/Genomics.pdf>`.


.. toctree::
   :numbered:
   :maxdepth: 2

   cli/index
   ngs-tools/index
   ngs-qc/index
   ngs-assembly/index
   ngs-annotation/index
   ngs-mapping/index
   ngs-taxonomic-investigation/index
   ngs-variantcalling/index
   general/quickref
   general/code
   general/downloads
   general/references
