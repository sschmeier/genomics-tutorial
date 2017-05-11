NGS - Genome annotation
=======================

Preface
-------

Overview
--------

The part of the workflow we will work on in this section can be viewed in :numref:`fig-workflow-anno`.

.. _fig-workflow-anno:
.. figure:: images/workflow.png

   The part of the workflow we will work on in this section marked in red.


Learning outcomes
-----------------

After studying this section of the tutorial you should be able to:

#. Explain how annotation completeness is assessed using orthologues
#. Use bioinformatics tools to perform gene prediction
#. Use genome-viewing software to graphically explore genome annotations and NGS data overlays 


Before we start
---------------

.. Attention:: The annotation process will take up to 90 minutes. Start it as soon
as possible.


Lets see how our directory structure looks so far:

.. code:: bash

          cd ~/analysis
          ls -1F

.. code:: bash

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



         
Genome annotation
-----------------
