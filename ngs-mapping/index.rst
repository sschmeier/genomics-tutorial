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

We are going to use a program called |bwa| fo map our reads to a genome.

It is simple to install and use.

.. code:: bash

          source activate ngs
          conda install samtools
          conda install bamtools
          conda install bwa

|bwa| overview
--------------

|bwa| is a short read aligner, that can take a reference genome and map single- or paired-end data to it.
It requires an indexing step in which one supplies the reference genome and |bwa| will create an index that in the subsequent steps will be used for aligning the reads to the reference genome.
Te general command structure of |bwa| looks like this:

.. code:: bash

   # bwa index help
   bwa index
          
   # indexing 
   bwa index reference-genome.fa

   # bwa mem help
   bwa mem
   
   # single-end mapping
   bwa mem reference-genome.fa reads.fq > aln-se.sam
   
   # paired-end mapping
   bwa mem reference-genome.fa read1.fq read2.fq > aln-pe.sam

   
Creating a reference index for mapping
--------------------------------------

.. todo::

   Create an |bwa| index for our reference genome assembly. Attention! Remember which file you need to submit to |bwa|.


.. hint::

   Should you not get it right, try these commands `here <../_static/code/mapping1.txt>`__.


Mapping reads in a paired-end manner
------------------------------------

Now that we have created our index, it is time to map the filtered and trimmed sequencing reads of our evolved line to the reference genome.

.. todo::
   
   Use the correct ``bwa mem`` command structure from above and map the reads of the evolved line to the reference genome.
   

.. hint::

   Should you not get it right, try these commands `here <../_static/code/mapping2.txt>`__.

   
The sam mapping file-format
---------------------------

|bwa| will produce a mapping file in sam-format. Have a look into the sam-file that was created by |bwa|. A quick overview of the sam-format can be found `here <http://bio-bwa.sourceforge.net/bwa.shtml#4>`__ and even more information can be found `here <http://samtools.github.io/hts-specs/SAMv1.pdf>`__. Briefly, first there are a lot of header lines. Then, for each read, that mapped to the reference, there is one line.

The columns of such a line in the mapping file are:

+-----+---------+-----------------------------------------------------------+
| Col |  Field	| Description                                               |
+=====+=========+===========================================================+
| 1   | QNAME   | Query (pair) NAME                                         |
+-----+---------+-----------------------------------------------------------+
| 2   | FLAG    | bitwise FLAG                                              |
+-----+---------+-----------------------------------------------------------+
| 3   | RNAME   | Reference sequence NAME                                   |
+-----+---------+-----------------------------------------------------------+
| 4   | POS     | 1-based leftmost POSition/coordinate of clipped sequence  |
+-----+---------+-----------------------------------------------------------+
| 5   | MAPQ    | MAPping Quality (Phred-scaled)                            |
+-----+---------+-----------------------------------------------------------+
| 6   | CIAGR   | extended CIGAR string                                     |
+-----+---------+-----------------------------------------------------------+
| 7   | MRNM    | Mate Reference sequence NaMe (‘=’ if same as RNAME)       |
+-----+---------+-----------------------------------------------------------+
| 8   | MPOS    | 1-based Mate POSition                                     |
+-----+---------+-----------------------------------------------------------+
| 9   |	ISIZE   | Inferred insert SIZE                                      |
+-----+---------+-----------------------------------------------------------+
| 10  |	SEQ     | query SEQuence on the same strand as the reference        |
+-----+---------+-----------------------------------------------------------+
| 11  |	QUAL    | query QUALity (ASCII-33 gives the Phred base quality)     |
+-----+---------+-----------------------------------------------------------+
| 12  |	OPT     | variable OPTional fields in the format TAG\:VTYPE\:VALUE  |
+-----+---------+-----------------------------------------------------------+

One line of a mapped read can be seen here:

.. code:: bash

    M02810:197:000000000-AV55U:1:1101:10000:11540   83      NODE_1_length_1419525_cov_15.3898       607378  60      151M    =       607100  -429    TATGGTATCACTTATGGTATCACTTATGGCTATCACTAATGGCTATCACTTATGGTATCACTTATGACTATCAGACGTTATTACTATCAGACGATAACTATCAGACTTTATTACTATCACTTTCATATTACCCACTATCATCCCTTCTTTA FHGHHHHHGGGHHHHHHHHHHHHHHHHHHGHHHHHHHHHHHGHHHHHGHHHHHHHHGDHHHHHHHHGHHHHGHHHGHHHHHHFHHHHGHHHHIHHHHHHHHHHHHHHHHHHHGHHHHHGHGHHHHHHHHEGGGGGGGGGFBCFFFFCCCCC NM:i:0  MD:Z:151        AS:i:151        XS:i:0

It basically defines, the read and the position in the reference genome where the read mapped and a quality of the map.

Sorting and compressing
-----------------------

We are going to use |samtools| to sort the sam-file and create a binary version for efficient storing of and access to the mapped reads. We are going to do the transformation into a bam-file (the binary version of a sam-file) and the sorting in one step:

.. rst-class:: sebcode

    samtools view -bS mappings/|fileevol|.sam | samtools sort -o mappings/|fileevol|.sorted.bam

.. attention::

   This step might take a few minutes to finish, depending on how big your mapping file is.
    

Mapping statistics
------------------

.. rst-class:: sebcode

    samtools depth mappings/|fileevol|.sorted.bam > mappings/|fileevol|.depth.txt

Unmapped reads
--------------

We could decide to use |kraken| like in section :ref:`taxonomic-investigation`
to classify all unmapped sequence reads and identify the species they are coming
from and test for contamination.
