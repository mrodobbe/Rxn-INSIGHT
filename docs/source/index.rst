.. Rxn-INSIGHT documentation master file, created by
   sphinx-quickstart on Thu Dec  5 15:53:43 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Rxn-INSIGHT documentation
=========================

Rxn-INSIGHT is an open-source algorithm designed to automate reaction classification, naming,
and the identification of functional groups, rings, and scaffolds from chemical entities.
It streamlines synthesis planning by analyzing reaction databases to recommend solvents, catalysts, and reagents
for novel reactions with high accuracy and efficiency. The tool replicates the reasoning of organic chemists,
achieving over 90% accuracy in classification and 95% in naming,
while providing suggestions for reaction conditions in under a second.

Rxn-INSIGHT can be installed from PyPI using pip (on python 3.10 or 3.11):

``pip install rxn_insight``

Please cite the Rxn-INSIGHT `paper <https://doi.org/10.1186/s13321-024-00834-z>`__, if you use it in your own work.

.. footbibliography::

.. toctree::
   :maxdepth: 2
   :caption: Quickstart

   installation
   getting-started
   advanced-features

.. toctree::
   :maxdepth: 2
   :caption: Guides

   common-use-cases
   visualization-guide

.. toctree::
   :maxdepth: 2
   :caption: API Documentation

   rxn_insight


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
