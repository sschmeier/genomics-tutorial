NGS - Quality control
=====================

Preface
-------

In this quality control section we will use our skill on the
command-line interface to deal with the task of investigating the quality and cleaning sequencing data.

There is an accompanying lectures for this tutorial:

-  `Next-generation sequencing and quality control: An introduction <https://dx.doi.org/10.6084/m9.figshare.2972320.v1>`__
   is available at `figshare <https://dx.doi.org/10.6084/m9.figshare.2972320.v1>`__
   (`https://dx.doi.org/10.6084/m9.figshare.2972320.v1 <https://dx.doi.org/10.6084/m9.figshare.2972320.v1>`__).

.. NOTE::

   You will encounter some **To-do** sections at times. Write the solutions and answers into a text-file.


Learning outcomes
-----------------

After studying this tutorial you should be able to:

#. Describe the steps involved in pre-processing/cleaning sequencing
   data.
#. Distinguish between a good and a bad sequencing run.
#. Compute, investigate and evaluate the quality of sequence data from a
   sequencing experiment.
   

The data
--------

First, we are going to download the data we will analyse. Open a shell/terminal.

.. code-block:: bash

   # create a directory you work in
   mkdir analysis

   # change into the directory
   cd analysis

   # download the data
   curl -O http://compbio.massey.ac.nz/data/203341/data.tar.gz

   # uncompress it
   tar -xvzf data.tar.gz


The data is from a paired-end sequencing run data (see :numref:`fig-pairedend`), thus we have two files, one
for each end of the read. 

.. _fig-pairedend:
.. figure:: images/pairedend.png

   Illustration of single-end (SE) versus paired-end (PE) sequencing.

If you need to refresh how Illumina paired-end sequencing works have a
look at the `Illumina
webpage <http://www.illumina.com/technology/next-generation-sequencing/paired-end-sequencing_assay.html>`__
and this `video <https://youtu.be/HMyCqWhwB8E>`__.

.. note::

   The data we are using is "almost" raw data coming from the machine. This data has been post-processed in two ways already. All sequences that were identified as belonging to the PhiX genome have been removed. This process requires some skills we will learn in later sections. Illumina adapters have been removed as well already! The process is explained below but we are not going to do it.


Investigate the data
~~~~~~~~~~~~~~~~~~~~

Make use of your newly developed skills on the command-line to
investigate the files in ``data`` folder.

.. admonition:: TODO

   #. Use the command-line to get some ideas about the file.
   #. What kind of files are we dealing with?
   #. How many sequence reads are in the file?


The QC process
--------------

There are a few steps one need to do when getting the raw sequencing data from the sequencing facility:

#. Remove PhiX sequences
#. Adapter trimming
#. Quality trimming of reads
#. Quality assessment
   

PhiX genome
-----------

`PhiX <https://en.wikipedia.org/wiki/Phi_X_174>`__ is a nontailed bacteriophage with a single-stranded DNA and a genome with 5386 nucleotides.
PhiX is used as a quality and calibration control for `sequencing runs <http://www.illumina.com/products/by-type/sequencing-kits/cluster-gen-sequencing-reagents/phix-control-v3.html>`__.
PhiX is often added at a low known concentration, spiked in the same lane along with the sample or used as a separate lane.
As the concentration of the genome is known, one can calibrate the instruments.
Thus, PhiX genomic sequences need to be removed before processing your data further as this constitutes a deliberate contamination [MUKHERJEE2015]_.
The steps involve mapping all reads to the "known" PhiX genome, and removing all of those sequence reads from the data.

However, your sequencing provider might not have used PhiX, thus you need to read the protocol carefully, or just do this step in any case.

.. note::

   We are not going to do this step here, as this has been already done. Please see the :doc:`../ngs-mapping/index` section on how to map reads against a genome.


Adapter trimming
----------------

The process of sequencing DNA via Illumina technology requires the addition of some adapters to the sequences.
These get sequenced as well and need to be removed as they are artificial and do not belong to the species we try to sequence.

.. note::

   The process of how to do this is explained here, however we are not going to do this as our sequences have been already adapter-trimmed.
   

Install a tool called `fastq-mcf <https://github.com/ExpressionAnalysis/ea-utils/blob/wiki/FastqMcf.md>`__  from the `ea-utils suite <https://expressionanalysis.github.io/ea-utils/>`__ of tools that is able to do this.

.. code-block:: bash

   # install
   conda install ea-utils


.. todo::

   Write this section.
   

Quality assessment of sequencing reads (|solexaqa|)
---------------------------------------------------

To assess the sequence read quality of the Illumina run we make use of a program called |solexaqa| [COX2010]_.
|solexaqa| was originally developed to work with Solexa data (since bought by Illumina), but long since working with Illumina data.
It produces nice graphics that intuitively show the quality of the sequences. it is also able to dynamically trim the bad quality ends off the reads.

From the webpage:

    "SolexaQA calculates sequence quality statistics and creates visual
    representations of data quality for second-generation sequencing
    data. Originally developed for the Illumina system (historically
    known as "Solexa"), SolexaQA now also supports Ion Torrent and 454
    data."

    
Install SolexaQA++
~~~~~~~~~~~~~~~~~~

Unfortunately, currently we cannot install |solexaqa| with |conda|.

.. code:: bash

    curl -O http://compbio.massey.ac.nz/data/203341/SolexaQA.tar.gz
   
    # uncompress the archive
    tar -xvzf SolexaQA.tar.gz
    
    # make the file executable
    chmod a+x SolexaQA/Linux_x64/SolexaQA++

    # copy program to root folder
    cp ./SolexaQA/Linux_x64/SolexaQA++ .
    
    # run the program
    ./SolexaQA++


.. note::

   Should the download fail, download manually from :doc:`../general/downloads`.

    
SolexaQA++ manual
~~~~~~~~~~~~~~~~~

|solexaqa| has three modes that can be run. Type:

.. code:: bash

     ./SolexaQA++
     
.. code:: bash

    SolexaQA++ v3.1.3
    Released under GNU General Public License version 3
    C++ version developed by Mauro Truglio (M.Truglio@massey.ac.nz)

    Usage: SolexaQA++ <command> [options]

    Command: analysis      quality analysis and graphs generation
             dynamictrim    trim reads using a chosen threshold
             lengthsort  sort reads by a chosen length

The three modes are: ``analysis``, ``dynamictrim``, and ``lengthsort``:

``analysis`` - the primary quality analysis and visualization tool.
Designed to run on unmodified FASTQ files obtained directly from
Illumina, Ion Torrent or 454 sequencers.

``dynamictrim`` - a read trimmer that individually crops each read to
its longest contiguous segment for which quality scores are greater than
a user-supplied quality cutoff.

``lengthsort`` - a program to separate high quality reads from low
quality reads. LengthSort assigns trimmed reads to paired-end, singleton
and discard files based on a user-defined length cutoff.


SolexaQA++ dynamic trimming
~~~~~~~~~~~~~~~~~~~~~~~~~~~

We will use |solexaqa| dynamic trim the reads, to chop of nucleotides witha a bad quality score.

.. admonition:: TODO

    #. Create a directory for the result-files --> **trimmed/**.
    #. Run |solexaqa| ``dynamictrim`` with the untrimmed data and a probability cutoff of 0.01., and submit result-directory **trimmed/**.
    #. Investigate the result-files in **trimmed/**, e.g. do the file-sizes change to the original files?
    #. |solexaqa| ``dynamictrim`` produces a graphical output. Explain what the graph shows. Find heklp on the |solexaqa| website.

.. hint::

   Should you not get 1 and/or 2 right, try these commands `here <../_static/code/qc-solexaqa1.txt>`__.

   
SolexaQA++ analysis on trimmed data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


.. admonition:: TODO

    #. Create a directory for the result-files --> **trimmed-solexaqa**.
    #. use |solexaqa| to do the quality assessment with the trimmed data-set.
    #. Compare your results to the examples of a particularly bad MiSeq run (:numref:`solexaqa_heatmap_bad` to :numref:`solexaqa_heatmap_bad`, taken from |solexaqa| website). Write down your observations.
    #. What elements in these example figures (:numref:`solexaqa_quality_bad` to :numref:`solexaqa_heatmap_bad`) indicate that the show a bad run? Write down your explanations.

.. hint::

   Should you not get 1 and/or 2 it right, try these commands `here <../_static/code/qc-solexaqa2.txt>`__ 


.. _solexaqa_quality_bad:
.. figure:: images/solexaqa_quality_bad.png

   SolexaQA++ example quality plot along reads of a bad MiSeq run

.. _solexaqa_hist_bad:
.. figure:: images/solexaqa_hist_bad.png

   SolexaQA++ example histogram plot of a bad MiSeq run.

.. _solexaqa_cumulative_bad:
.. figure:: images/solexaqa_cumulative_bad.png

   SolexaQA++ example cumulative plot of a bad MiSeq run.

.. _solexaqa_heatmap_bad:
.. figure:: images/solexaqa_heatmap_bad.png

   SolexaQA++ example quality heatmap of a bad MiSeq run.


   
|sickle| for dynamic trimming (alternative to |solexaqa|)
---------------------------------------------------------


Should the dynamic trimming not work with |solexaqa|, you can alternatively use |sickle|.

.. code:: bash

    source activate ngs
    conda install sickle-trim

Now we are going to run the program on our paired-end data:

.. rst-class:: sebcode

    # create a new directory
    mkdir trimmed
    
    # sickle parameters:
    sickle --help

    # as we are dealing with paired-end data you will be using "sickle pe"
    sickle pe --help

    # run sickle like so:
    sickle pe -g -t sanger -f data/|fileanc1|.fastq.gz -r data/|fileanc2|.fastq.gz -o trimmed/|fileanc1|.trimmed.fastq.gz -p trimmed/|fileanc2|.trimmed.fastq.gz 
  

.. hint::

   Should you be unable to run |sickle| or |solexaqa| at all to trim the data. You can download the trimmed dataset `here <http://compbio.massey.ac.nz/data/203341/trimmed.tar.gz>`__. Unarchive and uncompress the files with ``tar -xvzf trimmed.tar.gz``.


Quality assessment of sequencing reads (|fastqc|)
-------------------------------------------------

      
Installing |fastqc|
~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    source activate ngs   
    conda install fastqc

    # should now run the program
    fastqc --help
    

.. code:: bash


                FastQC - A high throughput sequence QC analysis tool

    SYNOPSIS

            fastqc seqfile1 seqfile2 .. seqfileN

        fastqc [-o output dir] [--(no)extract] [-f fastq|bam|sam]
               [-c contaminant file] seqfile1 .. seqfileN

    DESCRIPTION

        FastQC reads a set of sequence files and produces from each one a quality
        control report consisting of a number of different modules, each one of
        which will help to identify a different potential type of problem in your
        data.

        If no files to process are specified on the command line then the program
        will start as an interactive graphical application.  If files are provided
        on the command line then the program will run with no user interaction
        required.  In this mode it is suitable for inclusion into a standardised
        analysis pipeline.

        
FastQC manual
~~~~~~~~~~~~~

|fastqc| is a very simple program to run that provides similar and additional information to |solexaqa|.

From the webpage:

    "FastQC aims to provide a simple way to do some quality control
    checks on raw sequence data coming from high throughput sequencing
    pipelines. It provides a modular set of analyses which you can use
    to give a quick impression of whether your data has any problems of
    which you should be aware before doing any further analysis."

    
The basic command looks like:

.. code:: bash

    $ fastqc -o RESULT-DIR INPUT-FILE.[txt/fa/fq] ...

-  ``-o RESULT-DIR`` is the directory where the result files will be written
-  ``INPUT-FILE.[txt/fa/fq]`` is the sequence file to analyze, can be more than one file.

.. hint::

   The result will be a HTML page per input file that can be opened in a web-browser.


Run FastQC on the untrimmed and trimmed data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. admonition:: TODO

   #. Create a directory for the results --> **trimmed-fastqc**
   #. Run FastQC on all **trimmed** files.
   #. Visit the |fastqc| website and read about sequencing QC reports for good and bad Illumina sequencing runs.
   #. Compare your results to these examples (:numref:`fastqc-bad1` to :numref:`fastqc-bad3`) of a particularly bad run (taken from the |fastqc| website) and write down your observations with regards to your data.
   #. What elements in these example figures (:numref:`fastqc-bad1` to :numref:`fastqc-bad3`) indicate that the example is from a bad run?

.. hint::

   Should you not get it right, try these commands `here <../_static/code/qc-fastqc.txt>`__.

.. _fastqc-bad1:
.. figure:: images/fastqc_bad1.png

    Quality score across bases.

.. _fastqc-bad2:
.. figure:: images/fastqc_bad2.png
            
    Quality per tile.

.. _fastqc-bad3:
.. figure:: images/fastqc_bad3.png
            
    GC distribution over all sequences.


Web links
---------

- Lectures for this topic: `Next-generation sequencing and quality control: An introduction <https://dx.doi.org/10.6084/m9.figshare.2972320.v1>`__
- `Illumina paired-end sequencing <http://www.illumina.com/technology/next-generation-sequencing/paired-end-sequencing_assay.html>`__
- |fastqc|
- |solexaqa|

