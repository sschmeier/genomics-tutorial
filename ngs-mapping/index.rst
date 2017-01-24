.. _ngs-mapping:

NGS - Read mapping
==================

Preface
-------

In this section we will use our skill on the command-line interface to map our
reads from the evolved line to our ancestral reference genome.

There is an accompanying lecture for this tutorial:

-  `Genome Assembly: An Introduction <https://dx.doi.org/10.6084/m9.figshare.2972323.v1>`__ available at
   `figshare <https://dx.doi.org/10.6084/m9.figshare.2972323.v1>`__
   (https://dx.doi.org/10.6084/m9.figshare.2972323.v1).

.. NOTE::

   You will encounter some **To-do** sections at times. Write the solutions and answers into a text-file.   


Learning outcomes
-----------------

After studying this tutorial you should be able to:

#. Explain the process of sequence read mapping.
#. Use bioinformatics tools to map sequencing reads to a reference genome.
#. Filter mapped reads based on quality.


Before we start
---------------

Lets see how our directory structure looks so far:

.. code:: bash

          cd ~/analysis
          ls -1F

.. code:: bash
          
          annotation/
          assembly/
          data/
          kraken/
          SolexaQA/
          SolexaQA++
          trimmed/
          trimmed-fastqc/
          trimmed-solexaqa/
          

Mapping sequence reads to a reference genome
--------------------------------------------

We want to map the sequencing reads to the ancestral reference genome we created in the section :ref:`ngs-assembly`.
We are going to use the quality trimmed forward and backward DNA sequences of the evolved line and use a program called |bwa| to map the reads.

.. todo::
                
   #. Discuss briefly why we are using the ancestral genome as a reference genome as opposed to a genome for the evolved line.

      
Installing the software
~~~~~~~~~~~~~~~~~~~~~~~

We are going to use a program called |spades| fo assembling our genome.
In a recent evaluation of assembly software, |spades| was found to be a good
choice for fungal genomes [ABBAS2014]_.
It is also simple to install and use.

.. code:: bash

          source activate ngs
          conda install samtools
          conda install bamtools
          conda install bwa


Creating a reference index for mapping
--------------------------------------


Mapping reads in a paired-end manner
------------------------------------


The sam mapping file-format
---------------------------


Unmapped reads
--------------

We could decide to use |kraken| like in section :ref:`taxonomic-investigation`
to classify all unmapped sequence reads and identify the species they are coming
from and test for contamination.
