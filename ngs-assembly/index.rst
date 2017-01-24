.. _ngs-assembly:

NGS - Genome assembly
=====================

Preface
-------

In this section we will use our skill on the command-line interface to create a
genome assembly from sequencing data.

There is an accompanying lecture for this tutorial:

-  `Genome Assembly: An Introduction <https://dx.doi.org/10.6084/m9.figshare.2972323.v1>`__ available at
   `figshare <https://dx.doi.org/10.6084/m9.figshare.2972323.v1>`__
   (https://dx.doi.org/10.6084/m9.figshare.2972323.v1).

.. NOTE::

   You will encounter some **To-do** sections at times. Write the solutions and answers into a text-file.   


Learning outcomes
-----------------

After studying this tutorial you should be able to:

#. Compute and interpret a whole genome assembly.
#. Judge the quality of a genome assembly.


Before we start
---------------

Lets see how our directory structure looks so far:

.. code:: bash

          cd ~/analysis
          ls -1F

.. code:: bash
          
          data/
          kraken/
          SolexaQA/
          SolexaQA++
          trimmed/
          trimmed-fastqc/
          trimmed-solexaqa/
          

Creating a genome assembly
--------------------------

We want to create a genome assembly for our ancestor.
We are going to use the quality trimmed forward and backward DNA sequences and
use a program called |spades| to build a genome assembly.

.. todo::
                
   #. Discuss briefly why we are using the ancestral sequences to create a
      reference genome as opposed to the evolved line.

      
Installing the software
~~~~~~~~~~~~~~~~~~~~~~~

We are going to use a program called |spades| fo assembling our genome.
In a recent evaluation of assembly software, |spades| was found to be a good
choice for fungal genomes [ABBAS2014]_.
It is also simple to install and use.

.. code:: bash

          source activate ngs
          conda install spades

          
|spades| usage
~~~~~~~~~~~~~~

.. code:: bash

    # change to your analysis root folder
    cd ~/analysis
    
    # first create a output directory for the assemblies
    mkdir assembly
    
    # to get a help for spades and an overview of the parameter type:
    spades.py -h


The two files we need to submit to |spades| are two paired-end read files.

.. rst-class:: sebcode
               
    spades.py -o assembly/spades_default/ -1 trimmed/|fileanc1|.fastq.trimmed.gz -2 trimmed/|fileanc2|.fastq.trimmed.gz                   


.. todo::
   
   #. Run |spades| with default parameters on the ancestor
   #. Read in the |spades| manual about about assembling with 2x150bp reads
   #. Run |spades| a second time but use the options suggested at the |spades| manual `section 3.4 <http://spades.bioinf.spbau.ru/release3.9.1/manual.html#sec3.4>`__ for assembling 2x150bp paired-end reads (are fungi multicellular?). Use a different output directory ``assembly/spades_150`` for this run.

.. hint::

   Should you not get it right, try these commands `here <../_static/code/assembly1.txt>`__.

   
Assembly quality assessment
---------------------------

Assembly statistics
~~~~~~~~~~~~~~~~~~~

|quast| (QUality ASsesment Tool) [GUREVICH2013]_, evaluates genome assemblies by computing various metrics, including:

-  N50: length for which the collection of all contigs of that length or
   longer covers at least 50% of assembly length
-  NG50: where length of the reference genome is being covered
-  NA50 and NGA50: where aligned blocks instead of contigs are taken
-  missassemblies: misassembled and unaligned contigs or contigs bases
-  genes and operons covered

It is easy with |quast| to compare these measures among several assemblies.
The program can be used on their website (`http://quast.bioinf.spbau.ru/
<http://quast.bioinf.spbau.ru/>`__).

We can install it locally with:

.. code:: bash

          source activate ngs
          conda install quast

Run |quast| with both assembly scaffolds.fasta files to compare the results.

.. hint::

   Should you be unable to run |spades| on the data, you can download the assemblies `here <http://compbio.massey.ac.nz/data/203341/assembly.tar.gz>`__. Unarchive and uncompress the files with ``tar -xvzf assembly.tar.gz``.


.. rst-class:: sebcode

   quast -o assembly/quast assembly/spades_default/scaffolds.fasta assembly/spades_150/scaffolds.fasta
   

.. todo::
               
   #. Compare the results of |quast| with regards to the two different assemblies.
   #. Which one do you prefer and why?
   
      
Assemblathon
------------

.. todo::
                
   Now that you know the basics for assembling a genome and judging their quality, play with the |spades| parameters to create the best assembly possible. 
    
   
Further reading
---------------

Background on Genome Assemblies
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  How to apply de Bruijn graphs to genome assembly. [COMPEAU2011]_ 
-  Sequence assembly demystified. [NAGARAJAN2013]_ 

Evaluation of Genome Assembly Software
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- GAGE: A critical evaluation of genome assemblies and assembly algorithms. [SALZBERG2012]_ 
- Assessment of de novo assemblers for draft genomes: a case study with fungal genomes. [ABBAS2014]_




Web links
---------

- Lectures for this topic: `Genome Assembly: An Introduction <https://dx.doi.org/10.6084/m9.figshare.2972323.v1>`__
- |spades| 
- `Quast <http://quast.bioinf.spbau.ru/>`__
- `Bandage <https://rrwick.github.io/Bandage/>`__ (Bioinformatics Application for Navigating De novo Assembly Graphs Easily) is a program that visualizes a genome assembly as a graph [WICK2015]_.

