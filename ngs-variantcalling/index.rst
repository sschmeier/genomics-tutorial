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
          conda install tabix

          
Preprocessing
-------------

We first need to make an index of our reference genome as this is required by the SNP caller.
Given a scaffold/contig file in fasta-format, e.g. ``scaffolds.fasta`` which is
located ina  directory ``assembly/spades_final``, use |samtools| to do this:

.. code:: bash
          
          samtools faidx assembly/spades_final/scaffolds.fasta
   

Furthermore we need to pre-process our mapping files a bit further and creaee a bam-index file (``.bai``) for each o the bam-files, e.g.:

.. rst-class:: sebcode
               
          bamtools index -in mappings/|fileevol|.bam


Lets also create a new directory for the variants:

.. code:: bash

          mkdir variants

Calling variants
----------------
          
|samtools| mpileup
~~~~~~~~~~~~~~~~~~

We use the sorted bam-file that we produced in the mapping step before.

.. rst-class:: sebcode

          # We first pile up all the reads
          samtools mpileup -g -f assembly/spades_final/scaffolds.fasta mappings/|fileevol|.sorted.bam > variants/|fileevol|.mpileup.bcf
          # Now we call the variants
          bcftools view -c -v variants/|fileevol|.mpileup.bcf > variants/|fileevol|.mpileup.vcf

          
|Freebayes|
~~~~~~~~~~~

Now we can do some variant calling with another tool called |freebayes|.
Given a reference genome scaffold file in fasta-format, e.g. ``scaffolds.fasta`` and the index in ``.fai`` format and a mapping file (e.g. "|fileevol|.sorted.bam") and a mapping index, we can call |freebayes| like so:

.. rst-class:: sebcode

          # Now we call variants and pipe the results into a new file
          freebayes -f assembly/spades_final/scaffolds.fasta mappings/|fileevol|.sorted.bam > variants/|fileevol|.freebayes.vcf

          
This will result in a variants file "|fileevol|.freebayes.vcf".


Post-processing
---------------

Understanding the output files (.vcf)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Lets look at a vcf-file:

.. rst-class:: sebcode


          # first 10 lines, which are part of the header
          cat variants/|fileevol|.freebayes.vcf | head

          
.. code:: bash
          
          ##fileformat=VCFv4.1
          ##fileDate=20161122
          ##source=freeBayes v1.0.2-29-g41c1313
          ##reference=genome/scaffolds.fasta
          ##contig=<ID=NODE_1_length_1394677_cov_15.3771,length=1394677>
          ##contig=<ID=NODE_2_length_1051867_cov_15.4779,length=1051867>
          ##contig=<ID=NODE_3_length_950567_cov_15.4139,length=950567>
          ##contig=<ID=NODE_4_length_925223_cov_15.3905,length=925223>
          ##contig=<ID=NODE_5_length_916389_cov_15.4457,length=916389>
          ##contig=<ID=NODE_6_length_772252_cov_15.4454,length=772252>

Lets look at the variants:

.. rst-class:: sebcode
               
          # remove header lines and look at top 4 entires
          cat variants/|fileevol|.freebayes.vcf | egrep -v '##' | head -4

          
.. code:: bash
          
          #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  unknown
          NODE_1_length_1394677_cov_15.3771       137621  .       T       C       76.5197 .       AB=0.318182;ABP=9.32731;AC=1;AF=0.5;AN=2;AO=7;CIGAR=1X;DP=22;DPB=22;DPRA=0;EPP=18.2106;EPPR=4.31318;GTI=0;LEN=1;MEANALT=1;MQM=56.1429;MQMR=56.4;NS=1;NUMALT=1;ODDS=17.6193;PAIRED=1;PAIREDR=1;PAO=0;PQA=0;PQR=0;PRO=0;QA=268;QR=540;RO=15;RPL=0;RPP=18.2106;RPPR=6.62942;RPR=7;RUN=1;SAF=7;SAP=18.2106;SAR=0;SRF=12;SRP=14.7363;SRR=3;TYPE=snp       GT:DP:DPR:RO:QR:AO:QA:GL    0/1:22:22,7:15:540:7:268:-17.3644,0,-42.2185
          NODE_1_length_1394677_cov_15.3771       568696  .       G       A       1269.62 .       AB=0;ABP=0;AC=2;AF=1;AN=2;AO=38;CIGAR=1X;DP=38;DPB=38;DPRA=0;EPP=3.23888;EPPR=0;GTI=0;LEN=1;MEANALT=1;MQM=60;MQMR=0;NS=1;NUMALT=1;ODDS=57.2844;PAIRED=1;PAIREDR=0;PAO=0;PQA=0;PQR=0;PRO=0;QA=1438;QR=0;RO=0;RPL=20;RPP=3.23888;RPPR=0;RPR=18;RUN=1;SAF=20;SAP=3.23888;SAR=18;SRF=0;SRP=0;SRR=0;TYPE=snp      GT:DP:DPR:RO:QR:AO:QA:GL        1/1:38:38,38:0:0:38:1438:-129.701,-11.4391,0
          NODE_1_length_1394677_cov_15.3771       612771  .       T       C       60.7485 .       AB=0.3;ABP=9.95901;AC=1;AF=0.5;AN=2;AO=6;CIGAR=1X;DP=20;DPB=20;DPRA=0;EPP=4.45795;EPPR=8.59409;GTI=0;LEN=1;MEANALT=1;MQM=49.5;MQMR=54.3571;NS=1;NUMALT=1;ODDS=13.9879;PAIRED=1;PAIREDR=1;PAO=0;PQA=0;PQR=0;PRO=0;QA=223;QR=540;RO=14;RPL=6;RPP=16.0391;RPPR=33.4109;RPR=0;RUN=1;SAF=4;SAP=4.45795;SAR=2;SRF=4;SRP=8.59409;SRR=10;TYPE=snp    GT:DP:DPR:RO:QR:AO:QA:GL        0/1:20:20,6:14:540:6:223:-12.5734,0,-40.0605


Statistics and filter
~~~~~~~~~~~~~~~~~~~~~

Now we can use it to do some statistics and filter our variant calls.
          
.. rst-class:: sebcode
               
          # get statistics
          rtg vcfstats variants/|fileevol|.freebayes.vcf

          
.. rst-class:: sebcode
          
          # only keep entries with qual of min 30
          rtg vcffilter -q 30 -i variants/|fileevol|.freebayes.vcf -o variants/|fileevol|.freebayes-q30.vcf

          
.. rst-class:: sebcode
          
          # look at stats for filtered
          rtg vcfstats variants/|fileevol|.freebayes-q30.vcf
          

Finding variants of interest (VAI)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Things to consider when looking for VAI:

- The quality score of the variant call.
  
  * Do we call the variant with a higher then normal score?
    
- The mapping quality score.
  
  * How confident are we that the reads where mapped here correctly?
- The location of the SNP.
  
  * SNPs in larger contigs probably more interesting than in tiny contigs.
  * Does the SNP overlap a coding region in the genome annotation?
    
- The type of SNP.
 
          
Overlap variants with genes
~~~~~~~~~~~~~~~~~~~~~~~~~~~
