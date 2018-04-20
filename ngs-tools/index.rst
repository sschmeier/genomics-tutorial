.. _tool-installation:

Tool installation
=================

Install the conda package manager
---------------------------------

We will use the package/tool managing system |conda| to install some programs
that we will use during the course.
It is not installed by default, thus we need to install it first to be able to use it. 


.. code-block:: bash

    # download latest conda installer
    curl -O https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh

    # run the installer
    bash Miniconda3-latest-Linux-x86_64.sh
    
    # delete the installer after successful run
    rm Miniconda3-latest-Linux-x86_64.sh


.. Note::
   Should the conda installer download fail. Please find links to alternative locations on the
   :doc:`../downloads` page.

    
Update ``.bashrc`` and ``.zshrc`` config-files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Before we are able to use |conda| we need to tell our shell where it can find the program.
We add the right path to the |conda| installation to our shell config files:

.. code::
   
   echo 'export PATH="/home/manager/miniconda3/bin:$PATH"' >> ~/.bashrc
   echo 'export PATH="/home/manager/miniconda3/bin:$PATH"' >> ~/.zshrc


.. Attention::
   The above assumes that your username is "manager", which is the default on a Biolinux install.
   Replace "manager" with your actual username.
   Find out with ``whoami``.
   

So what is actually happening here? We are appending a line to a file (either ``.bashrc`` or ``.zshrc``).
If we are starting a new command-line shell, either file gets executed first (depending on which shell you are using, either bash or zsh shells).
What this line does, is to put permanently the directory ``~/miniconda3/bin`` first on your ``PATH`` variable.
The ``PATH`` variable contains directories in which our computer looks for installed programs, one directory after the other until the program you requested is found (or not, then it will complain).
Through the addition of the above line we make sure that the program ``conda`` can be found anytime we open a new shell.


Close shell/terminal, **re-open** new shell/terminal.
Now, we should be able to use the |conda| command:


.. code-block:: bash

    conda update conda


Installing conda channels to make tools available
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Different tools are packaged in what |conda| calls channels.
We need to add some channels to make the bioinformatics and genomics tools
available for installation:


.. code-block:: bash
    
    # Install some conda channels
    # A channel is where conda looks for packages
    conda config --add channels defaults
    conda config --add channels conda-forge
    conda config --add channels bioconda      

   
Create environments
-------------------

We create a |conda| environment for some tools This is useful to work **reproducible** as we can easily re-create the tool-set with the same version numbers later on.


.. code-block:: bash

    conda create -n ngs python=3
    # activate the environment
    conda activate ngs

    
So what is happening when you type ``conda activate ngs`` in a shell.
The ``PATH`` variable (mentioned above) gets temporarily manipulated and set to:


.. code-block:: bash
                
   $ conda activate ngs
   # Lets look at the content of the PATH variable
   (ngs) $ echo $PATH
   /home/manager/miniconda3/envs/ngs/bin:/home manager/miniconda3/bin:/usr/local/bin: ...


Now it will look first in your environment's bin directory but afterwards in the general conda bin (/home/manager/miniconda3/bin).
So basically everything you install generally with conda (without being in an environment) is also available to you but gets overshadowed if a similar program is in ``/home/manager/miniconda3/envs/ngs/bin`` and you are in the ``ngs`` environment.


Install software
----------------

To install software into the activated environment, one uses the command ``conda install``.

.. code-block:: bash
         
    # install more tools into the environment
    conda install package


.. note::
   To tell if you are in the correct conda environment, look at the command-prompt.
   Do you see the name of the environment in round brackets at the very beginning of the prompt, e.g. (ngs)?
   If not, activate the ``ngs`` environment with ``conda activate ngs`` before installing the tools.

    
                
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
    conda activate [name]

    # deavtivate env
    conda deactivate
