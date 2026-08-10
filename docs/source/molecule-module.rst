Molecule Module
================

The enhanced Molecule module extends Rxn-INSIGHT's capabilities with PubChem integration, similarity calculations, and additional analysis features. This allows for more comprehensive molecular analysis and information retrieval.

Features
--------

- **PubChem Integration:** Retrieve chemical names, descriptions, and identifiers
- **Similarity Calculation:** Compare molecules using Morgan fingerprints and Tanimoto similarity
- **Functional Group Analysis:** Identify chemical functional groups
- **Ring System Detection:** Extract ring structures from molecules
- **Reaction Search:** Find reactions that produce or involve the molecule

Installation
-------------

The enhanced Molecule module requires additional dependencies:

.. code-block:: bash

    pip install requests tqdm

Basic Usage
------------

Creating a Molecule and Analyzing Properties
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    import rxn_insight as ri

    # Create a molecule from SMILES
    mol = ri.Molecule("c1ccc(C(=O)O)cc1")  # Benzoic acid

    # Get basic properties
    print(f"SMILES: {mol.smiles}")
    print(f"InChI: {mol.inchi}")
    print(f"InChIKey: {mol.inchikey}")

    # Identify functional groups
    functional_groups = mol.get_functional_groups()
    print(f"Functional groups: {functional_groups}")

    # Find ring systems
    rings = mol.get_rings()
    print(f"Ring systems: {rings}")

    # Get molecular scaffold
    print(f"Scaffold: {mol.scaffold}")

Creating a Molecule from a Chemical Name
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``Compound`` is a ``Molecule`` built from a *name* rather than a SMILES string.
The name is resolved with an online service, after which the object supports
everything ``Molecule`` does:

.. code-block:: python

    import rxn_insight as ri

    compound = ri.Compound("benzene")

    print(compound.smiles)     # 'c1ccccc1'
    print(compound.scaffold)   # 'c1ccccc1'
    print(compound.inchikey)   # 'UHOVQNZJYSORNB-UHFFFAOYSA-N'

    # Everything Molecule offers is available
    print(compound.get_functional_groups())

By default the name is resolved through `OPSIN <https://opsin.ch.cam.ac.uk>`__,
which interprets systematic IUPAC names. Pass ``use_opsin=False`` to resolve
through PubChem instead, which is more forgiving of trivial and trade names:

.. code-block:: python

    aspirin = ri.Compound("aspirin", use_opsin=False)          # via PubChem
    ethanol = ri.Compound("ethanol", allow_pubchem=True)       # also fetch descriptions

.. note::
    Both routes require network access. A name that cannot be resolved raises
    ``ValueError`` on either route, so a ``Compound`` that is constructed
    successfully is always fully initialised:

    .. code-block:: python

        try:
            compound = ri.Compound(user_supplied_name)
        except ValueError as exc:
            print(f"Could not resolve that name: {exc}")

PubChem Integration
~~~~~~~~~~~~~~~~~~~~

Retrieve chemical information from PubChem:

.. code-block:: python

    # Create a molecule with PubChem information
    aspirin = ri.Molecule("CC(=O)OC1=CC=CC=C1C(=O)O", allow_pubchem=True)

    # Access PubChem information
    print(f"IUPAC Name: {aspirin.iupac_name}")
    print(f"Common Name: {aspirin.trivial_name}")
    print(f"PubChem CID: {aspirin.cid}")

    # Description may contain information about the molecule's uses and properties
    if aspirin.description:
        print(f"Description excerpt: {aspirin.description[:100]}...")

Calculating Molecular Similarity
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Compare molecules using fingerprint-based similarity:

.. code-block:: python

    # Create reference molecule
    benzene = ri.Molecule("c1ccccc1")

    # Compare with other molecules
    similarities = {
        "Toluene": benzene.calculate_similarity("Cc1ccccc1"),
        "Phenol": benzene.calculate_similarity("Oc1ccccc1"),
        "Naphthalene": benzene.calculate_similarity("c1ccc2ccccc2c1"),
        "Cyclohexane": benzene.calculate_similarity("C1CCCCC1")
    }

    # Print results
    for name, similarity in similarities.items():
        print(f"Similarity to {name}: {similarity:.3f}")

Searching for Reactions
~~~~~~~~~~~~~~~~~~~~~~~~~~

Find reactions that produce the molecule:

.. code-block:: python

    import pandas as pd
    from rxn_insight.molecule import Molecule

    # Load a reaction database
    df_rxns = pd.read_parquet("your_reaction_database.gzip")

    # Create a molecule
    biphenyl = Molecule("c1ccc(-c2ccccc2)cc1")

    # Find reactions that produce this molecule
    reactions = biphenyl.search_reactions(df_rxns)

    # Print the reactions
    if reactions is not None and len(reactions) > 0:
        print(f"Found {len(reactions)} reactions producing this molecule:")
        for i, (idx, row) in enumerate(reactions.head(3).iterrows()):
            print(f"\nReaction {i+1}:")
            print(f"SMILES: {row['REACTION']}")
            print(f"Class: {row['CLASS']}")
            print(f"Conditions: {row['SOLVENT']}, {row['CATALYST']}, {row['REAGENT']}")
    else:
        print("No reactions found for this molecule.")

Finding Similar Reactions by Scaffold
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Search for reactions with products having similar scaffolds:

.. code-block:: python

    # Find reactions with similar scaffolds
    similar_rxns = biphenyl.search_reactions_by_scaffold(
        df_rxns,
        threshold=0.6,
        max_return=10,
        fp="Morgan"
    )

    # Print the similar reactions
    if similar_rxns is not None and len(similar_rxns) > 0:
        print(f"Found {len(similar_rxns)} reactions with similar scaffolds:")
        for i, (idx, row) in enumerate(similar_rxns.head(3).iterrows()):
            print(f"\nSimilar reaction {i+1} (Similarity: {row['SIMILARITY']:.2f}):")
            print(f"Product: {row['REACTION'].split('>>')[1]}")
            print(f"Reaction: {row['REACTION']}")
            print(f"Class: {row['CLASS']}")

API Reference
---------------

.. py:class:: Molecule(smi, allow_pubchem=False)

   A class to handle and analyze molecular structures.

   :param str smi: SMILES string of the molecule
   :param bool allow_pubchem: Whether to fetch PubChem information (default: False)

   .. py:attribute:: mol

      RDKit molecule object

   .. py:attribute:: smiles

      SMILES representation of the molecule

   .. py:attribute:: inchi

      InChI identifier of the molecule

   .. py:attribute:: inchikey

      InChIKey identifier of the molecule

   .. py:attribute:: trivial_name

      Common name from PubChem (if available)

   .. py:attribute:: iupac_name

      IUPAC name from PubChem (if available)

   .. py:attribute:: description

      Description from PubChem (if available)

   .. py:attribute:: cid

      PubChem compound ID (if available)

   .. py:attribute:: functional_groups

      List of functional groups in the molecule

   .. py:attribute:: rings

      Ring structures in the molecule

   .. py:attribute:: scaffold

      Murcko scaffold of the molecule

   .. py:attribute:: maccs_fp

      MACCS fingerprint

   .. py:attribute:: morgan_fp

      Morgan fingerprint

   .. py:attribute:: reactions

      DataFrame of reactions involving this molecule

   .. py:method:: get_pubchem_information()

      Retrieves information about the molecule from PubChem's REST API.

   .. py:method:: calculate_similarity(smi)

      Calculates the chemical similarity between this molecule and another one.

      :param str smi: SMILES string of the molecule to compare with
      :return: Tanimoto similarity value (0-1, where 1 is identical)
      :rtype: float

   .. py:method:: search_reactions(df)

      Searches for reactions involving the molecule as a product.

      :param pandas.DataFrame df: The DataFrame to search for reactions
      :return: DataFrame containing matching reactions
      :rtype: pandas.DataFrame

   .. py:method:: search_reactions_by_scaffold(df, threshold=0.5, max_return=100, fp='MACCS')

      Searches for reactions based on scaffold similarity.

      :param pandas.DataFrame df: DataFrame containing reactions to search
      :param float threshold: Similarity threshold to apply (0-1)
      :param int max_return: Maximum number of reactions to return
      :param str fp: Type of fingerprint to use ('MACCS' or 'Morgan')
      :return: DataFrame of similar reactions, sorted by similarity
      :rtype: pandas.DataFrame

   .. py:method:: get_functional_groups(df=None)

      Identifies and returns the functional groups present in the molecule.

      :param pandas.DataFrame df: DataFrame containing functional group patterns
      :return: List of functional group names found in the molecule
      :rtype: list[str]

   .. py:method:: get_rings()

      Identifies and returns rings in the molecule.

      :return: List of ring SMILES strings found in the molecule
      :rtype: list[str]

.. py:class:: Compound(name, allow_pubchem=False, use_opsin=True)

   A :py:class:`Molecule` created from a chemical name instead of a SMILES
   string. The name is resolved with an online service, after which every
   ``Molecule`` attribute and method is available.

   :param str name: Chemical name to resolve, e.g. ``"benzene"``
   :param bool allow_pubchem: Whether to fetch additional PubChem information
      once the structure is known (default: False)
   :param bool use_opsin: Resolve through OPSIN; set to False to resolve
      through PubChem instead (default: True)

   .. py:attribute:: name

      The chemical name the compound was created from

   :raises ValueError: If the name cannot be resolved to a structure

   .. py:method:: read_from_opsin()

      Resolves the name with the OPSIN web service, which interprets
      systematic IUPAC names.

      :raises ValueError: If OPSIN returns no structure for the name

   .. py:method:: read_from_pubchem()

      Resolves the name with the PubChem REST API and records the CID when
      available. More forgiving of trivial and trade names than OPSIN.

      :raises ValueError: If PubChem returns no structure for the name

   .. note::
      Resolution requires network access. Connection failures surface as the
      underlying ``requests`` exception rather than ``ValueError``, so they stay
      distinguishable from a name that simply could not be resolved.
