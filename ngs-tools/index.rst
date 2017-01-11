.. _tool-installation:

NGS - Tool installation
=======================

Install the conda package manager
---------------------------------

We will use the package/tool managing system |conda| to install some programs
that we will use during the course. It is not installed by default, thus we need
to install it first to be able to use it. 

.. code-block:: bash

    # download latest conda installer
    curl -O https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh

    # run the installer
    bash Miniconda3-latest-MacOSX-x86_64.sh
    
    # delete the installer after successful run
    rm Miniconda3-latest-MacOSX-x86_64.sh

    # Install some conda channels
    # A channel is where conda looks for packages
    conda config --add channels conda-forge
    conda config --add channels defaults
    conda config --add channels r
    conda config --add channels bioconda

    
Close shell/terminal, **re-open** new shell/terminal.

.. code-block:: bash

    conda update conda

    
.. ATTENTION::
   Should the conda installer download fail. Please find links to alternative locations on the
   :doc:`../general/downloads` page.
   

Create environment
------------------

We create a |conda| environment for some tools This is useful to work
**reproducible** as we can easily re-create the tool-set with the same version
numbers later on.

.. code-block:: bash

    conda create -n ngs python=3 
    # activate the environment
    source activate ngs


Install software
----------------

To install software into the activated environment, one uses the command ``conda install``.

.. code-block:: bash
         
    # install more tools into the environment
    conda install package

    
General conda commands
----------------------

.. code-block:: bash

    # to search for packages
    conda search [package]
    
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
