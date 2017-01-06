NGS - Taxonomic investigation
=============================


Overview
--------

We want to investigate if there are sequences of other species in our collection
of sequences DNA pieces. We know most of them are from our species that we try
to assemble. However, lets investigate if we have possibly sequences from other
species sequenced maybe through contamination.

We will use the tool |kraken| to assign
taxonomic classifications to our sequence reads. Let us see if we can id some
sequences from other fungi or even bacteria.


Kraken
------

We will be using a tool called |kraken| [WOOD2014]_. This tool uses
k-mers to assign a taxonomic labels in form of |ncbitax| to the sequence (if
possible). The taxonomic label is assigned based on similar k-mer content of the
sequence in question to the k-mer content of reference genome sequence. The
result is a classification of the sequence in question to the most likely
taxonomic label. If the k-mer content is not similar to any genomic sequence in
the database used, it will not assign any taxonomic label.


Installation
------------

Use conda in the same fashion as before to install |kraken|:

.. code:: bash
          
   source activate ngs
   conda install kraken-all

   
Now we need to create or download a |kraken| database that can be used to assign
the taxonomic labels to sequences. We opt for downloading a pre-build database
from the |kraken| website:

.. code:: bash
          
   curl -O https://ccb.jhu.edu/software/kraken/dl/minikraken.tgz

   # alternatively we can use wget
   wget https://ccb.jhu.edu/software/kraken/dl/minikraken.tgz
   
   # once the download is finished, we need to extract the archive content
   # it will create a directory: "minikraken_20141208/"
   tar -xvzf minikraken.tgz


.. ATTENTION::
   Should the download fail. Please find links to alternative locations on the
   :doc:`../general/downloads` page.
   
   
   
Usage
-----

Now that we have installed |kraken| and downloaded and extracted the minikraken
database, we can attempt to investigate the sequences we got back from the
sequencing provider for other species as the one it should contain. We call the
|kraken| tool and specify the database and fasta-file with the sequences it
should use. 

.. rst-class:: sebcode
   
   kraken --db minikraken_20141208 |filebase|.fa > |filebase|.kraken

   
This may take a few minutes, depending on how many sequences we are going to
classify. The resulting content of the file |filebase|.kraken looks similar to
the following example:

.. include:: example-kraken.txt
   :literal:
   :end-line: 5

Here, the first column indicates if the sequence could be classified (C) or is
unclassified (U). The second column represents the original sequence id. The
third column is the |ncbitax| identifier that has been predicted for the sequence.

.. NOTE::
   The |kraken| manual can be accessed `here <http://ccb.jhu.edu/software/kraken/MANUAL.html>`__.


Investigate taxa
----------------

We can use the webpage `NCBI TaxIdentifier
<https://www.ncbi.nlm.nih.gov/Taxonomy/TaxIdentifier/tax_identifier.cgi>`__ to
quickly get the names to the taxonomy identifier. However, this is impractical
as we are dealing potentially with many sequences. |kraken| has some scripts
that help us understand our results better.

First, we generate a sample-wide report of all taxa found. This can be achieved
with the tool `kranken-report`.


.. rst-class:: sebcode 

   kraken-report --db minikraken_20141208 |filebase|.kraken | gzip > |filebase|.kraken.report.gz

   
The first few lines of thus a report are shown below.
   
.. include:: example-kraken-report.txt 
   :literal: 
   
   
To be able to compare the taxa content of one sample to another, we can create a report whose
structure is always the same, disregarding which taxa are found (obviously the
percentages and numbers will be different) for all samples. We will print out all
taxa (instead of only those found), ``--show-zeros`` option and sort them
according to taxa-ids (column 5), e.g. ``sort -n -k5``.


.. rst-class:: sebcode 

   kraken-report *--show-zeros* --db minikraken_20141208 |filebase|.kraken | **sort -n -k5** | gzip > |filebase|.kraken.report.sorted.gz

The report is not ordered according to taxa ids and contains all taxa in the
database, even if they have not been found in our sample and are thus zero:
   
.. include:: example-kraken-report2.txt 
   :literal: 
   
   
Second, for every sequence in our sample and its predicted taxonomic identifier,
we attach the taxonomic names.


.. rst-class:: sebcode 

   kraken-translate --mpa-format --db --db minikraken_20141208 |filebase|.kraken | gzip > |filebase|.kraken.names.gz


.. include:: example-kraken-names.txt 
   :literal: 
   :end-line: 5

   
References
----------
   
.. [WOOD2014] 
   Kraken: ultrafast metagenomic sequence classification using exact
   alignments. Derrick E WoodEmail author and Steven L Salzberg. Genome Biology,
   2014, 15:R46, `DOI: 10.1186/gb-2014-15-3-r46 <http://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-3-r46>`__.
