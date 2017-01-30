Coding solutions
================

QC
--

.. _code-qc1:

Code: |fastqc| 
~~~~~~~~~~~~~~

*Create directory:*

.. rst-class:: sebcode

   mkdir trimmed-fastqc


*Run FastQC:*

.. rst-class:: sebcode

   fastqc -o trimmed-fastqc trimmed/|fileanc1|.fastq.trimmed.gz trimmed/|fileanc2|.fastq.trimmed.gz trimmed/|fileevol1|.fastq.trimmed.gz trimmed/|fileevol2|.fastq.trimmed.gz


*Open html webpages:*

.. rst-class:: sebcode

   firefox trimmed-fastqc/\*.html


.. _code-qc2:

Code: |solexaqa| trimming
~~~~~~~~~~~~~~~~~~~~~~~~~

*Create directory for result-files:*

.. rst-class:: sebcode

   mkdir trimmed


*Run SolexaQA++:*

.. rst-class:: sebcode
               
   ./SolexaQA++ dynamictrim -p 0.01 -d trimmed/ data/|fileanc1|.fastq.gz
   
   ./SolexaQA++ dynamictrim -p 0.01 -d trimmed/ data/|fileanc2|.fastq.gz

   ./SolexaQA++ dynamictrim -p 0.01 -d trimmed/ data/|fileevol1|.fastq.gz

   ./SolexaQA++ dynamictrim -p 0.01 -d trimmed/ data/|fileevol2|.fastq.gz


.. _code-qc3:

Code: |solexaqa| qc
~~~~~~~~~~~~~~~~~~~

*Create directory for result-files:*

.. rst-class:: sebcode
               
   mkdir trimmed-solexaqa/

   
*Run SolexaQA++:*

.. rst-class:: sebcode
               
   ./SolexaQA++ analysis -d trimmed-solexaqa trimmed/|fileanc1|.fastq.trimmed.gz
   
   ./SolexaQA++ analysis -d trimmed-solexaqa trimmed/|fileanc2|.fastq.trimmed.gz

   ./SolexaQA++ analysis -d trimmed-solexaqa trimmed/|fileevol1|.fastq.trimmed.gz

   ./SolexaQA++ analysis -d trimmed-solexaqa trimmed/|fileevol2|.fastq.trimmed.gz


Assembly
--------

.. _code-assembly1:

Code: |spades| assembly
~~~~~~~~~~~~~~~~~~~~~~~

.. rst-class:: sebcode

   spades.py -o assembly/spades-150/ -k 21,33,55,77 --careful  -1 trimmed/|fileanc1|.fastq.trimmed.gz -2 trimmed/|fileanc2|.fastq.trimmed.gz


   
Mapping
-------

.. _code-bowtie1:

Code: |bowtie| indexing
~~~~~~~~~~~~~~~~~~~~~~~

*Build the index:*

.. rst-class:: sebcode

   bowtie2-build assembly/spades-final/scaffolds.fasta assembly/spades-final/scaffolds


.. _code-bowtie2:

Code: |bowtie| mapping
~~~~~~~~~~~~~~~~~~~~~~
   
*Map to the genome. Use a max fragemnt length of 1000 bp:*

.. rst-class:: sebcode

   bowtie2 -X 1000 -x assembly/spades-final/scaffolds -1 trimmed/|fileevol1|.fastq.trimmed.gz -2 trimmed/|fileevol2|.fastq.trimmed.gz -S mappings/|fileevol|.sam 

   
.. _code-bwa1: 

Code: |bwa| indexing 
~~~~~~~~~~~~~~~~~~~~

*Index the genome assembly:*

.. rst-class:: sebcode
               
   bwa index assembly/spades-final/scaffolds.fasta


.. _code-bwa2:

Code: |bwa| mapping 
~~~~~~~~~~~~~~~~~~~

*Run bwa mem:*

.. rst-class:: sebcode

   bwa mem assembly/spades-final/scaffolds.fasta trimmed/|fileevol1|.fastq.trimmed.gz trimmed/|fileevol2|.fastq.trimmed.gz > mappings/|fileevol|.sam 
