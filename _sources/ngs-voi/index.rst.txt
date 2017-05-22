NGS - Variants-of-interest
==========================

Preface
-------

In this section we will use our genome annotation of our reference and our genome variants in the evolved line to find variants that are interesting in terms of the observed biology.

.. NOTE::

    You will encounter some **To-do** sections at times. Write the solutions and answers into a text-file.   


Overview
--------

The part of the workflow we will work on in this section can be viewed in :numref:`fig-workflow-voi`.

.. _fig-workflow-voi:
.. figure:: images/workflow.png
    
    The part of the workflow we will work on in this section marked in red.
   
     
Learning outcomes
-----------------

After studying this section of the tutorial you should be able to:

#. Identify variants of interests.
#. Understand how the variants might affect the observed biology in the evolved line.


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
    mappings/
    phylogeny/
    SolexaQA/
    SolexaQA++
    trimmed/
    trimmed-fastqc/
    trimmed-solexaqa/
    variants/

  
General comments for identifying variants-of-interest
-----------------------------------------------------


Things to consider when looking for variants-of-interest:

- The quality score of the variant call.
  
  * Do we call the variant with a higher then normal score?
    
- The mapping quality score.
  
  * How confident are we that the reads were mapped at the position correctly?
    
- The location of the SNP.
  
  * SNPs in larger contigs are probably more interesting than in tiny contigs.
  * Does the SNP overlap a coding region in the genome annotation?
    
- The type of SNP.

  * substitutions vs. indels 


SnpEff
------

We will be using |snpeff| to annotate our identified variants. The tool will tell us on to which genes we should focus further analyses.


Installing software
~~~~~~~~~~~~~~~~~~~
  
Tools we are going to use in this section and how to install them if you not have done it yet.


.. code:: bash

    # activate the env
    source activate ngs
          
    # Install these tools into the conda environment
    # if not already installed
    conda install snpeff
    conda install genometools-genometools
  

Make a directory for the results (in your analysis directory) and change into
the directory:


.. code:: bash

    mkdir voi

    # change into the directory
    cd voi

         
Prepare SnpEff database
~~~~~~~~~~~~~~~~~~~~~~~

We need to create our own config-file for |snpeff|. Where is the ``snpEff.config``:


.. code:: bash

    find ~ -name snpEff.config
    /home/manager/miniconda3/envs/ngs/share/snpeff-4.3.1m-0/snpEff.config
    

This will give you the path to the ``snpEff.config``. It might be looking a bit different then the one shown here.

Make a local copy of the ``snpEff.config`` and then edit it with an editor of your choice:


.. code:: bash

    cp /home/manager/miniconda3/envs/ngs/share/snpeff-4.3.1m-0/snpEff.config .
    nano snpEff.config

          
Make sure the data directory path in the ``snpEff.config`` looks like this:


.. code:: bash

    data.dir = ./data/

          
There is a section with databases, which starts like this:


.. code:: bash

    #-------------------------------------------------------------------------------
    # Databases & Genomes
    #
    # One entry per genome version. 
    #
    # For genome version 'ZZZ' the entries look like
    #	ZZZ.genome              : Real name for ZZZ (e.g. 'Human')
    #	ZZZ.reference           : [Optional] Comma separated list of URL to site/s Where information for building ZZZ database was extracted.
    #	ZZZ.chrName.codonTable  : [Optional] Define codon table used for chromosome 'chrName' (Default: 'codon.Standard')
    #
    #-------------------------------------------------------------------------------


Add the following two lines in the database section underneath these header lines:


.. code:: bash

    # my yeast genome
    yeastanc.genome : WildYeastAnc

          
Now, we need to create a local data folder called ``./data/yeastanc``.


.. code:: bash

    # create folders
    mkdir -p ./data/yeastanc


Copy our genome assembly to the newly created data folder.
The name needs to be ``sequences.fa`` or ``yeastanc.fa``:


.. code:: bash
    
    cp ../assembly/spades-final/scaffolds.fasta ./data/yeastanc/sequences.fa
    gzip ./data/yeastanc/sequences.fa

    
Copy our genome annotation to the data folder.
The name needs to be ``genes.gff`` (or ``genes.gtf`` for gtf-files).


.. code:: bash

    cp ../annotation/your_new_fungus.gff ./data/yeastanc/genes.gff
    gzip ./data/yeastanc/genes.gff


Now we can build a new |snpeff| database:


.. code:: bash

    snpEff build -c snpEff.config -gff3 -v yeastanc


.. note::
   Should this fail, due to gff-format of the annotation, we can try to convert the gff to gtf:


.. code:: bash

    # using genometools
    gt gff3_to_gtf -gzip ../annotation/your_new_fungus.gff -o ./data/yeastanc/genes.gtf.gz


Now, we can use the gtf annotation top build the database:


.. code:: bash
          
    snpEff build -c snpEff.config -gtf22 -v yeastanc


SNP annotation
~~~~~~~~~~~~~~

Now we can use our new |snpeff| database to annotate some variants, e.g.:


.. code:: bash

    snpEff -c snpEff.config yeastanc ../variants/evolved-6.freebayes.filtered.vcf.gz > evolved-6.freebayes.filtered.anno.vcf


