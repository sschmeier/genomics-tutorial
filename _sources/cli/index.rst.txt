Command-line interface introduction
===================================

Preface
-------

This tutorial is based on a Linux/Unix *command-line*.
Using the *command-line* requires a Linux/Unix operating system.
The easiest way to try out a Linux system without actually installing it on your
computer is a `LiveCD <https://en.wikipedia.org/wiki/Live_CD>`__.
A LiveCD is a CD/DVD that you prepare (e.g. burn a Linux distribution on
it) and insert in your computer.
You would restart you computer and can run Linux from the CD/DVD without any installation requirements.
This is helpful for trying out a distribution of Linux, not for actual work.

Another route would be to use a virtual machine. A virtual computer that runs within your normal host system, e.g. Windows or MacOSX.
The software to create a virtual machine is free, e.g. `VirtualBox <https://www.virtualbox.org/>`__.

Common flavors of Linux ready for download are e.g. `Ubuntu <https://help.ubuntu.com/community/LiveCD>`__ or if you are
thinking of going the bioinformatics route,
`BioLinux <http://environmentalomics.org/bio-linux/>`__, which includes many pre-installed bioinformatics tools (this is also the distribution we will be using).

An accompanying lecture for this tutorial is available at `figshare <http://dx.doi.org/10.6084/m9.figshare.1506799>`__.

Learning outcomes
-----------------

#. Be able to operate comfortably the command-line.
#. Be able to navigate the unix directory structure on the command-line.
#. Be able to start command-line programs and getting help/information
   about programs.
#. Be able to investigate text files with command-line commands.
#. Be able to investigate the content of text-files on the command-line.
#. Be able to explain the concept of a unix pipe.

Introduction
------------

We will go through some introductory material explaining the shell syntax.
This is bash syntax but most of it will work on other shells (tcsh, zsh) as well.

What is a shell? Here I shamelessly quote `Wikipedia <https://goo.gl/g9x4tE>`__:

    "In computing, a shell is a user interface for access to an
    operating system's services. In general, operating system shells use
    either a command-line interface (**CLI**) or graphical user
    interface (GUI), depending on a computer's role and particular
    operation..."

    "**CLI** shells allow some operations to be performed faster in some
    situations, especially when a proper GUI has not been or cannot be
    created. However, they require the user to memorize all commands and
    their calling syntax, and also to learn the shell-specific scripting
    language, for example bash script."

The Ubuntu Linux desktop environment
------------------------------------

The default environment in Ubuntu is called Unity and is similar to other user interfaces found in Windows or MacOSX (:numref:`fig-desktop1`).

.. _fig-desktop1:
.. figure:: images/cli_Desktop1.png

   The BioLinux desktop environment Unity.


Some words regarding the Linux file-system
------------------------------------------

The directory structure in a Linux system is not much different from any other system you worked with, e.g. Windows, MacOSX.
It is essentially a tree structure (:numref:`fig-dir1`).

.. _fig-dir1:
.. figure::  images/cli_dir1.png

   Quick look at the directory tree structure on the command-line.

To navigate the file-system you can use a file-manager e.g. "Files" the default file manager in the Unity window manager used by BioLinux (:numref:`fig-dir2`).

.. _fig-dir2:
.. figure::  images/cli_dir2.png

   Quick look at the directory tree structure in the "Files" GUI.

However, on the command-line we navigate via commands and not via mouse clicks.
Why is it necessary to use the command-line in the first place?
Strictly speaking it is not, if you do not want to make use of programs on the command-line.
However, the power of the Linux system becomes only obvious once we learn to make use of the command-line, thus navigating the directory structure via commands is one of the **most important skills** for you to learn.

Open a terminal
---------------

Open a terminal window and you are are ready to go.
On your linux desktop find: **Application** --> **Accessories** --> **Terminal** (for Gnome environment) or type "Terminal" in the search box (:numref:`fig-shell-out0`).


.. _fig-shell-out0:        
.. figure:: images/cli_shell_out0.png

    Unity search bar.

:numref:`fig-shell-out1` shows an example of how a terminal window might look like (it is very easy to change its appearance).
You will se this window to execute the commands to work with files and biological data.
However, it is by no means restricted to "biological data", once you know how to handle the command-line many tasks based on files will be easily achieved using various programs available here.

.. _fig-shell-out1:
.. figure:: images/cli_shell_out1.png

    An example of a terminal window in Unity.

.. note:: The command-line prompt (e.g. ``$`` or ``>``) indicates that the shell/terminal is waiting for commands from us. These will be sent to the computer to execute. As long as you do not see the prompt, the computer is busy processing your request.

Proxy settings
--------------

You might encounter problems connecting to the internet.
So this is most likely the case if your university has restrictions in place to make the network more secure.
One of these measures to make a network more secure is a proxy.
However, we need the internet. Follow these steps to get connected.

Open system settings
^^^^^^^^^^^^^^^^^^^^

.. _fig-proxy1:
.. figure:: images/cli_proxy2.png

   Accessing the system settings.

Open network settings
^^^^^^^^^^^^^^^^^^^^^

.. _fig-proxy2:
.. figure:: images/cli_proxy2.png

   Accessing network settings.

Change proxy settings to automatic
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. _fig-proxy3:
.. figure:: images/cli_proxy3.png

   Changing proxy settings.

After changing the proxy settings to automatic you can open a web-browser and you should be asked for you network *username* and *password*.
After you typed those and hit *Enter*, you should be connected.

Let's get started
-----------------

Pre-2016 you would find my own introductory course in this place.
However, since 2016 we will be using the excellent material from the `Software Carpentry <http://software-carpentry.org>`__ Foundation.
I am a SWC affiliated volunteer instructor and we are teaching basic computer skills to scientists, with the goal of general computational up-skilling in the sciences.
The material is created in a collaborative manner and tested over and over in many workshops.

Please follow the link to the material at the `Software Carpentry <http://swcarpentry.github.io/shell-novice>`__ and have the webpage open.
Otherwise you will work in the terminal window aka the shell.

