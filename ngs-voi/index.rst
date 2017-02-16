NGS - Variants of interest identification
=========================================

Preface
-------

In this section we will use our genome annotation of the reference and our genome variants in the evolved line to find variants that are interesting in terms of the observed biology.

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
#. Understand how the variants might affect the observed biology.


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
        

Finding variants of interest (VOI)
----------------------------------


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
