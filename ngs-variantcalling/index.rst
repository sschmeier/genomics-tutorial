NGS - Variant calling
=====================

Preface
-------

In this section we will use our genome assembly based on the ancestor and call
genetic variants in the evolved line.

There is an accompanying lecture for this tutorial:

- `SNPs - GWAS - eQTLs introduction <http://dx.doi.org/10.6084/m9.figshare.1515026>`__

.. NOTE::

   You will encounter some **To-do** sections at times. Write the solutions and answers into a text-file.   

   
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

  
Learning outcomes
-----------------

After studying this tutorial you should be able to:

#. Use tools to call variants based on a reference genome.
#. Identify variants of interests.
#. Understand how the variants might affect the observed biology.

   
Installing necessary software
-----------------------------
  
Tools we are going to use in this section and how to intall them if you not have done it yet.

.. code:: bash

          # activate the env
          source activate ngs
          
          # Install these tools into the conda environment
          # if not already installed
          conda install samtools
          conda install bamtools
          conda install bcftools
          conda install freebayes
          conda install rtg-tools
          conda install bedtools

          
Preprocessing
-------------

We first need to make an index of our reference genome as this is required by the SNP caller.
Given a scaffold/contig file in fasta-format, e.g. ``scaffolds.fasta`` which is located in the directory ``assembly/spades_final``, use |samtools| to do this:


.. code:: bash
          
          samtools faidx assembly/spades-final/scaffolds.fasta
   

Furthermore we need to pre-process our mapping files a bit further and create a bam-index file (``.bai``) for the bam-file we want to work with:


.. rst-class:: sebcode
               
          bamtools index -in mappings/|fileevol|.sorted.concordant.q20.bam


Lets also create a new directory for the variants:


.. code:: bash

          mkdir variants

          
Calling variants
----------------

SAMtools mpileup
~~~~~~~~~~~~~~~~

We use the sorted filtered bam-file that we produced in the mapping step before.

.. rst-class:: sebcode

   # We first pile up all the reads and then call variants
   samtools mpileup -u -g -f assembly/spades-final/scaffolds.fasta mappings/|fileevol|.sorted.concordant.q20.bam | bcftools call -v -m -O z -o variants/|fileevol|.mpileup.vcf.gz
   
|samtools| mpileup parameter:

- ``-u``: uncompressed output
- ``-g``: generate genotype likelihoods in BCF format
- ``-f FILE``: faidx indexed reference sequence file
  
|bcftools| view parameter:

- ``-v``: output variant sites only
- ``-m``: alternative model for multiallelic and rare-variant calling
- ``-o``: output file-name
- ``-O z``: output type: 'z' compressed VCF

  
Freebayes
~~~~~~~~~

As an alternative we can do some variant calling with another tool called |freebayes|.
Given a reference genome scaffold file in fasta-format, e.g. ``scaffolds.fasta`` and the index in ``.fai`` format and a mapping file (.bam file) and a mapping index (.bai file), we can call variants with |freebayes| like so:

.. rst-class:: sebcode

   # Now we call variants and pipe the results into a new file
   freebayes -f assembly/spades-final/scaffolds.fasta mappings/|fileevol|.sorted.concordant.q20.bam | gzip > variants/|fileevol|.freebayes.vcf.gz

         
Post-processing
---------------

Understanding the output files (.vcf)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Lets look at a vcf-file:

.. rst-class:: sebcode

   # first 10 lines, which are part of the header
   zcat variants/|fileevol|.mpileup.vcf.gz | head

          
.. code:: bash
   
   ##fileformat=VCFv4.2
   ##FILTER=<ID=PASS,Description="All filters passed">
   ##samtoolsVersion=1.3.1+htslib-1.3.1
   ##samtoolsCommand=samtools mpileup -g -f assembly/spades-final/scaffolds.fasta -o variants/evolved-6.mpileup.bcf mappings/evolved-6.sorted.concordant.q20.bam
   ##reference=file://assembly/spades-final/scaffolds.fasta
   ##contig=<ID=NODE_1_length_1419525_cov_15.3898,length=1419525>
   ##contig=<ID=NODE_2_length_1254443_cov_15.4779,length=1254443>
   ##contig=<ID=NODE_3_length_972329_cov_15.3966,length=972329>
   ##contig=<ID=NODE_4_length_951685_cov_15.4231,length=951685>
   ##contig=<ID=NODE_5_length_925222_cov_15.39,length=925222>
   ##contig=<ID=NODE_6_length_916533_cov_15.4426,length=916533>

Lets look at the variants:

.. rst-class:: sebcode
               
   # remove header lines and look at top 4 entires
   zcat variants/|fileevol|.mpileup.vcf.gz | egrep -v '##' | head -4

          
.. code:: bash
          
   #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  mappings/evolved-6.sorted.concordant.q20.bam
   NODE_1_length_1419525_cov_15.3898       24721   .       T       C       164     .       DP=12;VDB=0.205941;SGB=-0.680642;MQ0F=0;AC=2;AN=2;DP4=0,0,12,0;MQ=40     GT:PL   1/1:191,36,0
   NODE_1_length_1419525_cov_15.3898       157033  .       AAGAGAGAGAGAGAGAGAGAGAGA        AAGAGAGAGAGAGAGAGAGAGA  39.3328  .       INDEL;IDV=6;IMF=0.146341;DP=41;VDB=0.0813946;SGB=-0.616816;MQSB=1;MQ0F=0;ICB=1;HOB=0.5;AC=1;AN=2;DP4=13,17,3,3;MQ=42     GT:PL   0/1:75,0,255
   NODE_1_length_1419525_cov_15.3898       162469  .       T       C       19.609  .       DP=16;VDB=0.045681;SGB=-0.511536;RPB=0.032027;MQB=0.832553;BQB=0.130524;MQ0F=0;ICB=1;HOB=0.5;AC=1;AN=2;DP4=13,0,3,0;MQ=39        GT:PL   0/1:54,0,155


The fields in a vcf-file are described in he table (:numref:`table-vcf`) below:

.. _table-vcf:
.. table:: The vcf-file format fields.

   +-----+-----------+--------------------------------------------------------------------------------------+
   | Col | Field     | Description                                                                          |
   +=====+===========+======================================================================================+
   | 1   | CHROM     | Chromosome name                                                                      |
   +-----+-----------+--------------------------------------------------------------------------------------+
   | 2   | POS       | 1-based position. For an indel, this is the position preceding the indel.            |
   +-----+-----------+--------------------------------------------------------------------------------------+
   | 3   | ID        | Variant identifier. Usually the dbSNP rsID.                                          |
   +-----+-----------+--------------------------------------------------------------------------------------+
   | 4   | REF       | Reference sequence at POS involved in the variant. For a SNP, it is a single base.   |
   +-----+-----------+--------------------------------------------------------------------------------------+
   | 5   | ALT       | Comma delimited list of alternative seuqence(s).                                     |
   +-----+-----------+--------------------------------------------------------------------------------------+
   | 6   | QUAL      | Phred-scaled probability of all samples being homozygous reference.                  |
   +-----+-----------+--------------------------------------------------------------------------------------+
   | 7   | FILTER    | Semicolon delimited list of filters that the variant fails to pass.                  |
   +-----+-----------+--------------------------------------------------------------------------------------+
   | 8   | INFO      | Semicolon delimited list of variant information.                                     |
   +-----+-----------+--------------------------------------------------------------------------------------+
   | 9   | FORMAT    | Colon delimited list of the format of individual genotypes in the following fields.  |
   +-----+-----------+--------------------------------------------------------------------------------------+ 
   | 10+ | Sample(s) | Individual genotype information defined by FORMAT.                                   |
   +-----+-----------+--------------------------------------------------------------------------------------+


          
Statistics and filter
~~~~~~~~~~~~~~~~~~~~~

Now we can use it to do some statistics and filter our variant calls.

First, to prepare out vcf-file for querying we need to index it with ``tabix``:

.. rst-class:: sebcode

   tabix -p vcf variants/|fileevol|.mpileup.vcf.gz


- ``-p vcf``: input format 


We can get some quick stats with ``rtg vcfstats``:


.. rst-class:: sebcode
               
   rtg vcfstats variants/|fileevol|.mpileup.vcf.gz

   
Example output from ``rtg vcfstats``:


.. code::

   Location                     : variants/evolved-6.mpileup.vcf.gz
   Failed Filters               : 0
   Passed Filters               : 516
   SNPs                         : 399
   MNPs                         : 0
   Insertions                   : 104
   Deletions                    : 13
   Indels                       : 0
   Same as reference            : 0
   SNP Transitions/Transversions: 1.87 (286/153)
   Total Het/Hom ratio          : 3.20 (393/123)
   SNP Het/Hom ratio            : 8.98 (359/40)
   MNP Het/Hom ratio            : - (0/0)
   Insertion Het/Hom ratio      : 0.30 (24/80)
   Deletion Het/Hom ratio       : 3.33 (10/3)
   Indel Het/Hom ratio          : - (0/0)
   Insertion/Deletion ratio     : 8.00 (104/13)
   Indel/SNP+MNP ratio          : 0.29 (117/399)
   

   
However, we can also run |bcftools| to extract more detailed statistics about our variant calls:
   

.. rst-class:: sebcode
               
   bcftools stats -F assembly/spades-final/scaffolds.fasta -s - variants/|fileevol|.mpileup.vcf.gz > variants/|fileevol|.mpileup.vcf.gz.stats


- ``-s -``: list of samples for sample stats, "-" to include all samples
- ``-F FILE``: faidx indexed reference sequence file to determine INDEL context

  
Now we take the stats and make some plots (e.g. :numref:`fig-vcfstats`) which are particular of interest if having multiple samples, as one can easily compare them. However, we are only working with one here:


.. rst-class:: sebcode
   
   mkdir variants/plots
   plot-vcfstats -p variants/plots/ variants/|fileevol|.vcf.gz.stats

   
- ``-p``: The output files prefix, add a slash to create new directory.
   

.. _fig-vcfstats:
.. figure:: images/vcfstats.png
            
    Example of ``plot-vcfstats`` output.


Next, we filter out low quality reads.
We only include variants that have quality > 30.


.. rst-class:: sebcode

   # use rtg vcfffilter
   rtg vcffilter -q 30 -i variants/|fileevol|.mpileup.vcf.gz -o variants/|fileevol|.mpileup.q30.vcf.gz


- ``-i FILE``: input file
- ``-o FILE``: output file
- ``-q FLOAT``: minimal allowed quality in output.
  
   
or use |bcftools|:

.. rst-class:: sebcode

   # or use bcftools
   bcftools filter -O z -o variants/|fileevol|.mpileup.q30.vcf.gz -i'%QUAL>=30' variants/|fileevol|.mpileup.vcf.gz
   # bcftools filter does not index output, so we need to do it again
   tabix -p vcf variants/|fileevol|.mpileup.q30.vcf.gz
      

- ``-i'%QUAL>=30'``: we only include variants that have been called with quality >= 30.


Quick stats for the filtered variants:
  
.. rst-class:: sebcode 
          
   # look at stats for filtered 
   rtg vcfstats variants/|fileevol|.mpileup.q30.vcf.gz
  
  
.. todo::
    
   Look at the statistics. One ratio that is mentioned in the statistics is transition transversion ratio (*ts/tv*). Explain what this ratio is and why the observed ratio makes sense. 


   
Finding variants of interest (VAI)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Things to consider when looking for VAI:

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

   SEB: Write this section.
