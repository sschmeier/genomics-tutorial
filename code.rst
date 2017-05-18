Coding solutions
================

QC
--

.. _code-qc1:

Code: FastQC 
~~~~~~~~~~~~

*Create directory:*

.. code:: bash

   mkdir trimmed-fastqc


*Run FastQC:*

.. code:: bash

   fastqc -o trimmed-fastqc trimmed/ancestor-R1.fastq.trimmed.gz trimmed/ancestor-R2.fastq.trimmed.gz trimmed/evolved-6-R1.fastq.trimmed.gz trimmed/evolved-6-R2.fastq.trimmed.gz


*Open html webpages:*

.. code:: bash

   firefox trimmed-fastqc/\*.html


.. _code-qc2:

Code: SolexaQA++ trimming
~~~~~~~~~~~~~~~~~~~~~~~~~

*Create directory for result-files:*

.. code:: bash

   mkdir trimmed


*Run SolexaQA++:*

.. code:: bash
               
   ./SolexaQA++ dynamictrim -p 0.01 -d trimmed/ data/ancestor-R1.fastq.gz
   
   ./SolexaQA++ dynamictrim -p 0.01 -d trimmed/ data/ancestor-R2.fastq.gz

   ./SolexaQA++ dynamictrim -p 0.01 -d trimmed/ data/evolved-6-R1.fastq.gz

   ./SolexaQA++ dynamictrim -p 0.01 -d trimmed/ data/evolved-6-R2.fastq.gz


.. _code-qc3:

Code: SolexaQA++ qc
~~~~~~~~~~~~~~~~~~~

*Create directory for result-files:*

.. code:: bash
               
   mkdir trimmed-solexaqa/

   
*Run SolexaQA++:*

.. code:: bash
               
   ./SolexaQA++ analysis -d trimmed-solexaqa trimmed/ancestor-R1.fastq.trimmed.gz
   
   ./SolexaQA++ analysis -d trimmed-solexaqa trimmed/ancestor-R2.fastq.trimmed.gz

   ./SolexaQA++ analysis -d trimmed-solexaqa trimmed/evolved-6-R1.fastq.trimmed.gz

   ./SolexaQA++ analysis -d trimmed-solexaqa trimmed/evolved-6-R2.fastq.trimmed.gz


Assembly
--------

.. _code-assembly1:

Code: SPAdes assembly (trimmed data)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: bash 

   spades.py -o assembly/spades-150/ -k 21,33,55,77 --careful -1 trimmed/ancestor-R1.fastq.trimmed.gz -2 trimmed/ancestor-R2.fastq.trimmed.gz 


.. _code-assembly2:
   
Code: SPAdes assembly (original data)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: bash 

   spades.py -o assembly/spades-original/ -k 21,33,55,77 --careful -1 data/ancestor-R1.fastq.gz -2 data/ancestor-R2.fastq.gz 


   
Mapping
-------

.. _code-bowtie1:

Code: Bowtie2 indexing
~~~~~~~~~~~~~~~~~~~~~~

*Build the index:*

.. code:: bash

   bowtie2-build assembly/spades-final/scaffolds.fasta assembly/spades-final/scaffolds


.. _code-bowtie2:

Code: Bowtie2 mapping
~~~~~~~~~~~~~~~~~~~~~~
   
*Map to the genome. Use a max fragemnt length of 1000 bp:*

.. code:: bash

   bowtie2 -X 1000 -x assembly/spades-final/scaffolds -1 trimmed/evolved-6-R1.fastq.trimmed.gz -2 trimmed/evolved-6-R2.fastq.trimmed.gz -S mappings/evolved-6.sam 

   
.. _code-bwa1: 

Code: BWA indexing 
~~~~~~~~~~~~~~~~~~~~

*Index the genome assembly:*

.. code:: bash
               
   bwa index assembly/spades-final/scaffolds.fasta


.. _code-bwa2:

Code: BWA mapping 
~~~~~~~~~~~~~~~~~~~

*Run bwa mem:*

.. code:: bash

   bwa mem assembly/spades-final/scaffolds.fasta trimmed/evolved-6-R1.fastq.trimmed.gz trimmed/evolved-6-R2.fastq.trimmed.gz > mappings/evolved-6.sam 
