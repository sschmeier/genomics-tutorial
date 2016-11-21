NGS - Tool installation
=======================

Install the conda package manager
---------------------------------

We will use the package/tool managing system
`conda <http://conda.pydata.org/miniconda.html>`__ to install some
programs that we will use during the course. It is not installed by
default, thus we need to install it first to be able to use it.

.. code:: bash

    # download latest conda installer
    curl -O https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
    # run the installer
    bash Miniconda3-latest-MacOSX-x86_64.sh
    # delete the installer after successful run
    rm Miniconda3-latest-MacOSX-x86_64.sh
    # copy config file to home dir
    cp data/.condarc ~

Close shell/terminal, **re-open** new shell/terminal.

.. code:: bash

    conda update conda

Create environment
------------------

We create a `conda <http://conda.pydata.org/miniconda.html>`__
environment for some tools This is useful to work **reproducible** as we
can easily re-create the tool-set with the same version numbers later
on.

.. code:: bash

    conda create -n ngs fastqc
    # activate the environment
    source activate ngs


Install software
----------------

.. code:: bash
          
    # install more tools into the environment
    conda install spades     # assembler
    conda install quast      # QA assemblies
    conda install bwa        # mapper
    conda install samtools   # Working with mapping files
    conda install bamtools   # Working with mapping files
    conda install freebayes  # SNP caller
    conda install bcftools   # working with SNP files
    conda install igvtools   # for visualisation

General conda commands
----------------------

.. code:: bash

    # To update all packages
    conda update --all --yes

    # List all packages installed
    conda list [-n env]

    # conda list environments
    conda env list

    # create new env
    conda create -n [name] package [package] ...

    # activate env
    source activate [name]

    # deavtivate env
    source deactivate
