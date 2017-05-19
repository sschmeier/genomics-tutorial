NGS - Variants-of-interest
==========================

.. warning::
   **THIS PART OF THE TUTORIAL IS CURRENTLY UNDER ACTIVE DEVELOPMENT SO EXPECT CONTENT
   TO CHANGE**

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
          SolexaQA/
          SolexaQA++
          trimmed/
          trimmed-fastqc/
          trimmed-solexaqa/
          variants/

   
Installing necessary software
-----------------------------
  
Tools we are going to use in this section and how to intall them if you not have done it yet.

.. code:: bash

          # activate the env
          source activate ngs
          
          # Install these tools into the conda environment
          # if not already installed
        

Identification of variants-of-interest (VOI)
--------------------------------------------


Things to consider when looking for VOI:

- The quality score of the variant call.
  
  * Do we call the variant with a higher then normal score?
    
- The mapping quality score.
  
  * How confident are we that the reads were mapped at the position correctly?
    
- The location of the SNP.
  
  * SNPs in larger contigs are probably more interesting than in tiny contigs.
  * Does the SNP overlap a coding region in the genome annotation?
    
- The type of SNP.

  * substitutions vs. indels 

    
Overlap variants with genes
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. todo::

   OLIN: Write this section.

   
Visualise variants on the reference genome
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. todo::

   OLIN: Write this section.
GET DATA
--------

.. code:: bash

    cp ~/projects_current/203341/assembly/spades-default/scaffolds.fasta .
    cp ~/projects_current/203341/variants/evolved-6.mpileup.q30.vcf.gz .

ANNOTATE SCAFFOLD GENES/ORFS WITH SNAP
--------------------------------------

I need a annotation to try snpEFF, SNAP was just quick. Get HMM for
yeast.

.. code:: bash

    curl -O https://raw.githubusercontent.com/hyphaltip/fungi-gene-prediction-params/master/params/SNAP/saccharomyces_cerevisiae_S288C.hmm
    curl -O https://raw.githubusercontent.com/hyphaltip/fungi-gene-prediction-params/master/params/SNAP/saccharomyces_cerevisiae_rm11-1a_1.hmm
    snap saccharomyces_cerevisiae_S288C.hmm scaffolds.fasta -gff | gzip > annotation.gff.gz

This produces gff format. However, ideally I want gtf for building
snpEFF database and thats why I think it broke the build process below.

PREPARE SNPEFF DATABASE
-----------------------

Where is snpEff.config:

.. code:: bash

    find ~ -name snpEff.config
    /Users/sschmeie/miniconda2/envs/snpeff/share/snpeff-4.3.1m-0/snpEff.config
    /Users/sschmeie/miniconda2/pkgs/snpeff-4.3.1m-0/share/snpeff-4.3.1m-0/snpEff.config

The first one is the riht one as I have an env "snpeff" where I
installed snpEff into.

Make a local copy of the config file. Edit the config file. There is a
section with databases:

``bash  cp /Users/sschmeie/miniconda2/envs/snpeff/share/snpeff-4.3.1m-0/snpEff.config .  emacs snpEff.config``

Make sure data directory looks like this:

data.dir = ./data/

Add the following two lines in the database section:

my genome yeast
===============

yeast1.genome : Yeast

Now we need to create a local data folder called './data/yeast1'.

.. code:: bash

    # create folders
    mkdir -p ./data/yeast1

    # Copy genome to newly created folder, name needs to be sequences.fa or yeast1.fa
    cp scaffolds.fasta ./data/yeast1/sequences.fa
    gzip ./data/yeast1/sequences.fa

    # copy annotation to folder, name needs to be genes.gff.gz (or genes.gtf.gz for gtf-files)
    cp annotation.gff.gz ./data/yeast1/genes.gff.gz

    # now build db
    snpEff build -c snpEff.config -gff3 -v yeast1

    # or for gtf
    snpEff build -c snpEff.config -gtf22 -v yeast1

I tried the process with a gff file created above and it failed. If you
have gtf it might work. If you sent me yours + your genome I try.

Use snpEff for annotation with local config file:

.. code:: bash

    snpEff -c snpEff.config ...

