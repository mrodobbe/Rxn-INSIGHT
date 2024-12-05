.. _installation:

Installation
============

Rxn-INSIGHT can be installed from PyPI via pip or directly from the GitHub repository.

Option 1: Installing from PyPI using pip
----------------------------------------
.. code-block::

    conda create -n rxn-insight python=3.11
    conda activate rxn-insight
    pip install rxn-insight

.. note::
    Currently only python 3.10 and python 3.11 are supported.

Option 2: Installing from source using pip
------------------------------------------
.. code-block::

    conda create -n rxn-insight python=3.11
    conda activate rxn-insight
    git clone https://github.com/mrodobbe/Rxn-INSIGHT.git
    cd Rxn-INSIGHT
    pip install -e .

.. note::
    You can also use this option to install additional optional dependencies for development purposes,
    which are required to run the tests and build the docs by running ``pip install -e ".[test,doc]"``.

