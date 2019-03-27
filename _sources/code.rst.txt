.. _ngs-code:

Coding solutions
================


QC
--

.. _code-sickle:

Code: Sickle
~~~~~~~~~~~~

.. code:: bash

   # run sickle like this on the ancestor:
   sickle pe -g -t sanger -f data/ancestor-R1.fastq.gz -r data/ancestor-R2.fastq.gz -o trimmed/ancestor-R1.trimmed.fastq.gz -p trimmed/ancestor-R2.trimmed.fastq.gz -s trimmed/ancestor-singles.fastq.gz

   # run the evolved samples through sickle
   sickle pe -g -t sanger -f data/evolved-6-R1.fastq.gz -r data/evolved-6-R2.fastq.gz -o trimmed/evolved-6-R1.trimmed.fastq.gz -p trimmed/evolved-6-R2.trimmed.fastq.gz -s trimmed/evolved-6-singles.fastq.gz


.. _code-qc1:

Code: FastQC
~~~~~~~~~~~~

*Create directory:*

.. code:: bash

   mkdir trimmed-fastqc


*Run FastQC:*

.. code:: bash

   fastqc -o trimmed-fastqc trimmed/ancestor-R1.trimmed.fastq.gz trimmed/ancestor-R2.trimmed.fastq.gz trimmed/evolved-6-R1.trimmed.fastq.gz trimmed/evolved-6-R2.trimmed.fastq.gz


*Open html webpages:*

.. code:: bash

   firefox trimmed-fastqc/\*.html



Assembly
--------

.. _code-assembly1:

Code: SPAdes assembly (trimmed data)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: bash

   spades.py -o assembly/spades-150/ -k 21,33,55,77 --careful -1 trimmed/ancestor-R1.trimmed.fastq.gz -2 trimmed/ancestor-R2.trimmed.fastq.gz


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

   bowtie2-build assembly/spades_final/scaffolds.fasta assembly/spades_final/scaffolds


.. _code-bowtie2:

Code: Bowtie2 mapping
~~~~~~~~~~~~~~~~~~~~~~

*Map to the genome. Use a max fragemnt length of 1000 bp:*

.. code:: bash

   bowtie2 -X 1000 -x assembly/spades_final/scaffolds -1 trimmed/evolved-6-R1.trimmed.fsatq.gz -2 trimmed/evolved-6-R2.trimmed.fastq.gz -S mappings/evolved-6.sam


.. _code-bwa1:

Code: BWA indexing
~~~~~~~~~~~~~~~~~~~~

*Index the genome assembly:*

.. code:: bash

   bwa index assembly/spades_final/scaffolds.fasta


.. _code-bwa2:

Code: BWA mapping
~~~~~~~~~~~~~~~~~~~

*Run bwa mem:*

.. code:: bash

   # trimmed data
   bwa mem assembly/spades_final/scaffolds.fasta trimmed/evolved-6-R1.trimmed.fastq.gz trimmed/evolved-6-R2.trimmed.fastq.gz > mappings/evolved-6.sam

   # raw data
   bwa mem assembly/spades_final/scaffolds.fasta data/evolved-6-R1.fastq.gz data/evolved-6-R2.fastq.gz > mappings/evolved-6.raw.sam
