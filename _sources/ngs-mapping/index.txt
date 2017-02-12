.. _ngs-mapping:

NGS - Read mapping
==================

Preface
-------

In this section we will use our skill on the command-line interface to map our
reads from the evolved line to our ancestral reference genome.

There is an accompanying lecture for this tutorial:

-  `Genome Assembly: An Introduction <https://dx.doi.org/10.6084/m9.figshare.2972323.v1>`__.

.. NOTE::

   You will encounter some **To-do** sections at times. Write the solutions and answers into a text-file.   

   
Overview
--------

.. _fig-workflow-map:
.. figure:: images/workflow.png

   The part of the workflow we will work on in this section.
   

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
          conda install bedtools
          conda install bowtie2
          conda install bwa

          
Bowtie2
-------

Overview
~~~~~~~~

|bowtie| is a short read aligner, that can take a reference genome and map single- or paired-end data to it.
It requires an indexing step in which one supplies the reference genome and |bowtie| will create an index that in the subsequent steps will be used for aligning the reads to the reference genome.
The general command structure of the |bowtie| tools we are going to use are shown below:


.. code:: bash

   # bowtie2 help
   bowtie2-build
          
   # indexing 
   bowtie2-build genome.fasta PATH_TO_INDEX_PREFIX

   # paired-end mapping
   bowtie2 -X 1000 -x PATH_TO_INDEX_PREFIX -1 read1.fq.gz -2 read2.fq.gz -S aln-pe.sam


- ``-X``: Adjust the maximum fragment size (length of paired-end alignments + insert size) to 1000bp. This might be useful if you do not know the exact insert size of your data. The |bowtie| default is set to 500 which is `often considered too short <http://lab.loman.net/2013/05/02/use-x-with-bowtie2-to-set-minimum-and-maximum-insert-sizes-for-nextera-libraries/>`__.
  

Creating a reference index for mapping
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. todo::

   Create an |bowtie| index for our reference genome assembly. Attention! Remember which file you need to submit to |bowtie|.


.. hint::

   Should you not get it right, try the commands in :ref:`code-bowtie1`.



Mapping reads in a paired-end manner
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Now that we have created our index, it is time to map the filtered and trimmed sequencing reads of our evolved line to the reference genome.

.. todo::
   
   Use the correct ``bowtie2`` command structure from above and map the reads of the evolved line to the reference genome.
   

.. hint::

   Should you not get it right, try the commands in :ref:`code-bowtie2`.

          
BWA
---

.. Attention::

   If the mapping did not succeed with |bowtie|. We can use the aligner |bwa| explained in this section. If the mapping with |bowtie| did work, you can jump this section.


Overview
~~~~~~~~

|bwa| is a short read aligner, that can take a reference genome and map single- or paired-end data to it.
It requires an indexing step in which one supplies the reference genome and |bwa| will create an index that in the subsequent steps will be used for aligning the reads to the reference genome.
The general command structure of the |bwa| tools we are going to use are shown below:

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
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. todo::

   Create an |bwa| index for our reference genome assembly. Attention! Remember which file you need to submit to |bwa|.


.. hint::

   Should you not get it right, try the commands in :ref:`code-bwa1`.


Mapping reads in a paired-end manner
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Now that we have created our index, it is time to map the filtered and trimmed sequencing reads of our evolved line to the reference genome.

.. todo::
   
   Use the correct ``bwa mem`` command structure from above and map the reads of the evolved line to the reference genome.
   

.. hint::

   Should you not get it right, try the commands in :ref:`code-bwa2`.

   
The sam mapping file-format
---------------------------

|bwa| will produce a mapping file in sam-format. Have a look into the sam-file that was created by |bwa|.
A quick overview of the sam-format can be found `here <http://bio-bwa.sourceforge.net/bwa.shtml#4>`__ and even more information can be found `here <http://samtools.github.io/hts-specs/SAMv1.pdf>`__.
Briefly, first there are a lot of header lines. Then, for each read, that mapped to the reference, there is one line.

The columns of such a line in the mapping file are described in :numref:`table-sam`.

.. _table-sam:
.. table:: The sam-file format fields.

   +-----+---------+-----------------------------------------------------------+
   | Col |  Field  | Description                                               |
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
   | 9   | ISIZE   | Inferred insert SIZE                                      |
   +-----+---------+-----------------------------------------------------------+
   | 10  | SEQ     | query SEQuence on the same strand as the reference        |
   +-----+---------+-----------------------------------------------------------+
   | 11  | QUAL    | query QUALity (ASCII-33 gives the Phred base quality)     |
   +-----+---------+-----------------------------------------------------------+
   | 12  | OPT     | variable OPTional fields in the format TAG\:VTYPE\:VALUE  |
   +-----+---------+-----------------------------------------------------------+

One line of a mapped read can be seen here:

.. code:: bash

    M02810:197:000000000-AV55U:1:1101:10000:11540   83      NODE_1_length_1419525_cov_15.3898       607378  60      151M    =       607100  -429    TATGGTATCACTTATGGTATCACTTATGGCTATCACTAATGGCTATCACTTATGGTATCACTTATGACTATCAGACGTTATTACTATCAGACGATAACTATCAGACTTTATTACTATCACTTTCATATTACCCACTATCATCCCTTCTTTA FHGHHHHHGGGHHHHHHHHHHHHHHHHHHGHHHHHHHHHHHGHHHHHGHHHHHHHHGDHHHHHHHHGHHHHGHHHGHHHHHHFHHHHGHHHHIHHHHHHHHHHHHHHHHHHHGHHHHHGHGHHHHHHHHEGGGGGGGGGFBCFFFFCCCCC NM:i:0  MD:Z:151        AS:i:151        XS:i:0

It basically defines, the read and the position in the reference genome where the read mapped and a quality of the map.


Mapping post-processing
-----------------------

Fix mates and compress
~~~~~~~~~~~~~~~~~~~~~~

Because aligners can sometimes leave unusual `SAM flag <http://bio-bwa.sourceforge.net/bwa.shtml#4>`__ information on SAM records, it is helpful when working with many tools to first clean up read pairing information and flags with |samtools|.
We are going to produce also compressed bam output for efficient storing of and access to the mapped reads.


.. rst-class:: sebcode
               
   samtools fixmate -O bam |fileevol|.sam |fileevol|.fixmate.bam

   
- ``-O bam``: specifies that we want compressed bam output


.. attention:: 

   The step of sam to bam-file conversion might take a few minutes to finish, depending on how big your mapping file is. 


We will be using the `SAM flag <http://bio-bwa.sourceforge.net/bwa.shtml#4>`__ information later below to extract specific alignments. 

.. hint::

   A very useful tools to explain flags can be found `here <http://broadinstitute.github.io/picard/explain-flags.html>`__.

      
Once we have bam-file, we can also delete the original sam-file as it requires too much space.
   
  
.. rst-class:: sebcode

   rm mappings/|fileevol|.sam


Sorting
~~~~~~~

We are going to use |samtools| again to sort the bam-file into coordinate order:


.. rst-class:: sebcode

    # convert to bam file and sort
    samtools sort -O bam -o |fileevol|.sorted.bam |fileevol|.fixmate.bam
    

- ``-o``: specifies the name of the output file.
- ``-O bam``: specifies that the output will be bam-format
    

Mapping statistics
------------------

Stats with SAMtools
~~~~~~~~~~~~~~~~~~~

Lets get an mapping overview:

.. rst-class:: sebcode

    samtools flagstat mappings/|fileevol|.sorted.bam

    
.. todo::

   Look at the mapping statistics and understand `their meaning
   <https://www.biostars.org/p/12475/>`__. Discuss your results.
   Explain why we may find mapped reads that have their mate mapped to a different chromosome/contig?
   Can they be used for something?
         
   
For the sorted bam-file we can get read depth for at all positions of the reference genome, e.g. how many reads are overlapping the genomic position.


.. rst-class:: sebcode

    samtools depth mappings/|fileevol|.sorted.bam | gzip > mappings/|fileevol|.depth.txt.gz


.. todo::

   Extract the depth values for contig 20 and load the data into R, calculate some statistics of our scaffold.

   
.. rst-class:: sebcode
   
   zcat mappings/evolved-6.depth.txt.gz | egrep '^NODE_20_' | gzip >  mappings/NODE_20.depth.txt.gz

   
Now we quickly use some |R| to make a coverage plot for contig NODE20.
Open a |R| shell by typing ``R`` on the command-line of the shell.
   
.. code:: R

   x <- read.table('mappings/NODE_20.depth.txt.gz', sep='\t', header=FALSE,  strip.white=TRUE)

   # Look at the beginning of x
   head(x)

   # calculate average depth
   mean(x[,3])
   # std dev
   sqrt(var(x[,3]))
   
   # mark areas that have a coverage below 20 in red
   plot(x[,2], x[,3], col = ifelse(x[,3] < 20,'red','black'), pch=19, xlab='postion', ylab='coverage')

   # to save a plot
   png('mappings/covNODE20.png', width = 1200, height = 500)
   plot(x[,2], x[,3], col = ifelse(x[,3] < 20,'red','black'), pch=19, xlab='postion', ylab='coverage')
   dev.off()


The result plot will be looking similar to the one in :numref:`coverage`

.. _coverage:
.. figure:: images/covNODE20.png

   A example coverage plot for a contig with highlighted in red regions with a coverage below 20 reads.
   
   
.. todo::

   Look at the created plot. Explain why it makes sense that you find relatively bad coverage at the beginning and the end of the contig.


Stats with QualiMap
~~~~~~~~~~~~~~~~~~~

For a more in depth analysis of the mappings, one can use |qualimap|.

|qualimap| examines sequencing alignment data in SAM/BAM files according to the features of the mapped reads and provides an overall view of the data that helps to the detect biases in the sequencing and/or mapping of the data and eases decision-making for further analysis.

Installation:

.. code::

   conda install qualimap
   
   
Sub-selecting reads
-------------------

It is important to remember that the mapping commands we used above, without additional parameter to sub-select specific alignments (e.g. for |bowtie| there are options like ``--no-mixed``, which suppresses unpaired alignments for paired reads or ``--no-discordant``, which suppresses discordant alignments for paired reads, etc.), is going to output all reads, including unmapped reads, multi-mapping reads, unpaired reads, discordant read pairs, etc. in one file. We can sub-select from the output reads we want to analyse further using |samtools|.


.. todo::

   Explain what concordant and discordant read pairs are? Look at the |bowtie| manual.
   

Concordant reads
~~~~~~~~~~~~~~~~

Here, we select the reads **we will be using for subsequent analyses**.
Frist off, we select reads with a mapping quality of at least 20.
Furthermore, we select read-pair that have been mapped in a correct manner (same chromosome/contig, correct orientation to each other).


.. rst-class:: sebcode
               
   samtools view -h -b -q 20 -f 2 mappings/|fileevol|.sorted.bam > mappings/|fileevol|.sorted.concordant.q20.bam


- ``-h``: Include the sam header
- ``-b``: Output will be bam-format
- ``-q 20``: Only extract reads with mapping quality >= 20
- ``-f 2``: Only extract correctly paired reads. ``-f`` extracts alignments with the specified `SAM flag <http://bio-bwa.sourceforge.net/bwa.shtml#4>`__ set.


.. attention::

   The resulting file of this step will be used in the next section for calling variants.


Unmapped reads
~~~~~~~~~~~~~~

We could decide to use |kraken| like in section :ref:`taxonomic-investigation` to classify all unmapped sequence reads and identify the species they are coming from and test for contamination.

Lets see how we can get the unmapped portion of the reads from the bam-file:


.. rst-class:: sebcode
               
    samtools view -b -f 4 mappings/|fileevol|.sorted.bam > mappings/|fileevol|.sorted.unmapped.bam
    
    # count them
    samtools view -c mappings/|fileevol|.sorted.unmapped.bam
    
    
- ``-b``: indicates that the output is BAM.
- ``-f INT``: only include reads with this `SAM flag <http://bio-bwa.sourceforge.net/bwa.shtml#4>`__ set. You can also use the command ``samtools flags`` to get an overview of the flags. 
- ``-c``: count the reads


Lets extract the fastq sequence of the unmapped reads for read1 and read2.


.. rst-class:: sebcode

    bamToFastq -i |fileevol|.sorted.unmapped.bam -fq mappings/|fileevol|.sorted.unmapped.R1.fastq -fq2  mappings/|fileevol|.sorted.unmapped.R2.fastq


