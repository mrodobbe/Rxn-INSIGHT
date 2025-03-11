Advanced Features and API Reference
===================================

This guide covers advanced features of the Rxn-INSIGHT package and
provides a concise API reference for the core classes.

Reaction Class
--------------

The ``Reaction`` class is the central component for analyzing chemical
reactions.

Key Attributes
~~~~~~~~~~~~~~

- ``reaction``: SMILES representation of the reaction
- ``reactants``: SMILES string of reactants
- ``products``: SMILES string of products
- ``mapped_reaction``: Reaction with atom mappings
- ``reaction_class``: Classification of the reaction
- ``name``: Name of the reaction
- ``scaffold``: Molecular scaffold of the product
- ``byproducts``: Tuple of byproducts from the reaction
- ``template``: Extracted reaction template

Important Methods
~~~~~~~~~~~~~~~~~

``get_reaction_info()``
^^^^^^^^^^^^^^^^^^^^^^^

Returns a comprehensive dictionary with reaction details: - Reaction
class and name - Functional groups in reactants and products - Ring
systems - Byproducts - Scaffold information - Atom mapping information

``find_neighbors(df, fp='MACCS', concatenate=True, max_return=100, threshold=0.3, broaden=False, full_search=False)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Finds similar reactions in a database: - ``df``: Pandas DataFrame
containing reaction data - ``fp``: Fingerprint type (‘MACCS’ or
‘Morgan’) - ``concatenate``: Whether to concatenate reactant and product
fingerprints - ``max_return``: Maximum number of results to return -
``threshold``: Similarity threshold (0-1) - ``broaden``: Use broader
search criteria - ``full_search``: Perform a full database search
(slower)

``suggest_conditions(df)``
^^^^^^^^^^^^^^^^^^^^^^^^^^

Suggests optimal conditions based on similar reactions: - ``df``: Pandas
DataFrame containing reaction data - Returns: A dictionary with
suggested solvent, catalyst, and reagent

``get_class()``
^^^^^^^^^^^^^^^

Determines and returns the reaction class.

``get_name()``
^^^^^^^^^^^^^^

Determines and returns the reaction name.

``get_byproducts()``
^^^^^^^^^^^^^^^^^^^^

Calculates and returns likely byproducts.

``get_scaffold()``
^^^^^^^^^^^^^^^^^^

Extracts and returns the molecular scaffold.

``get_rings_in_reactants()``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Identifies ring structures in reactants.

``get_rings_in_products()``
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Identifies ring structures in products.

Molecule Class
--------------

The ``Molecule`` class handles operations related to individual
molecules.

.. _key-attributes-1:

Key Attributes
~~~~~~~~~~~~~~

- ``mol``: RDKit molecule object
- ``smiles``: SMILES representation
- ``inchi``: InChI identifier
- ``inchikey``: InChIKey identifier
- ``scaffold``: Murcko scaffold of the molecule
- ``maccs_fp``: MACCS fingerprint
- ``morgan_fp``: Morgan fingerprint

.. _important-methods-1:

Important Methods
~~~~~~~~~~~~~~~~~

``get_functional_groups(df=None)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Identifies functional groups in the molecule.

``get_rings()``
^^^^^^^^^^^^^^^

Extracts ring structures from the molecule.

``search_reactions(df)``
^^^^^^^^^^^^^^^^^^^^^^^^

Finds reactions in the database where this molecule is a product.

``search_reactions_by_scaffold(df, threshold=0.5, max_return=100, fp='MACCS')``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Finds reactions with similar product scaffolds.

Database Class
--------------

The ``Database`` class manages collections of reactions.

Key Methods
~~~~~~~~~~~

``create_database_from_df(df, reaction_column, solvent_column='SOLVENT', reagent_column='REAGENT', catalyst_column='CATALYST', yield_column='YIELD', ref_column='REF')``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Creates a reaction database from a DataFrame: - ``df``: Input DataFrame
with reaction data - ``reaction_column``: Column containing reaction
SMILES - Other parameters: Specify column names for conditions

``create_database_from_csv(fname, reaction_column, ...)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Creates a database from a CSV file.

``save_to_parquet(fname)``
^^^^^^^^^^^^^^^^^^^^^^^^^^

Saves the database to a parquet file.

``get_class_distribution()``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Returns the distribution of reaction classes.

``get_name_distribution()``
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Returns the distribution of reaction names.

Utility Functions
-----------------

The ``utils`` module contains various helper functions:

Reaction Handling
~~~~~~~~~~~~~~~~~

- ``get_atom_mapping(rxn, rxn_mapper=None)``: Maps atoms in a reaction
- ``get_reaction_template(reaction, radius_reactants=2, radius_products=2)``:
  Extracts a reaction template
- ``sanitize_mapped_reaction(rxn)``: Cleans up a mapped reaction
- ``remove_atom_mapping(rxn, smarts=False)``: Removes atom mapping

Fingerprinting and Similarity
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- ``get_fp(rxn, fp='MACCS', concatenate=True)``: Gets a fingerprint for
  a reaction
- ``get_similarity(v1, v2, metric='jaccard')``: Calculates similarity
  between fingerprints
- ``maccs_fp(mol)``: Gets MACCS fingerprint for a molecule
- ``morgan_fp(mol)``: Gets Morgan fingerprint for a molecule

Scaffold Analysis
~~~~~~~~~~~~~~~~~

- ``get_scaffold(mol)``: Gets the Murcko scaffold of a molecule
- ``get_ring_systems(mol, include_spiro=False)``: Identifies ring
  systems

Ranking Functions
~~~~~~~~~~~~~~~~~

- ``get_solvent_ranking(df)``: Ranks solvents by frequency
- ``get_catalyst_ranking(df)``: Ranks catalysts by frequency
- ``get_reagent_ranking(df)``: Ranks reagents by frequency

Advanced Usage Examples
-----------------------

Custom Reaction Classification
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

   from rxn_insight.reaction import Reaction
   from rxn_insight.classification import ReactionClassifier

   # Create a reaction
   reaction_smiles = "CC(=O)OC1=CC=CC=C1>>OC1=CC=CC=C1.CC(=O)O"

   # Access the classifier directly for advanced analysis
   rxn = Reaction(reaction_smiles)
   classifier = rxn.classifier

   # Directly check classification properties
   print(f"Is functional group interconversion: {classifier.is_fgi()}")
   print(f"Is deprotection: {classifier.is_deprotection()}")
   print(f"Is protection: {classifier.is_protection()}")
   print(f"Is oxidation: {classifier.is_oxidation()}")
   print(f"Is reduction: {classifier.is_reduction()}")
   print(f"Is C-C coupling: {classifier.is_cc_coupling()}")

Working with Atom Mappings
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

   from rxn_insight.reaction import Reaction
   from rxnmapper import RXNMapper

   # Initialize RXNMapper
   rxn_mapper = RXNMapper()

   # Map a reaction
   rxn_smiles = "OB(O)c1ccccc1.Brc1ccccc1>>c1ccc(-c2ccccc2)cc1"
   mapped_rxn = rxn_mapper.get_attention_guided_atom_maps([rxn_smiles])[0]["mapped_rxn"]

   # Create a Reaction with the mapping
   rxn = Reaction(mapped_rxn, keep_mapping=True)

   # Get the reaction center
   reaction_center = rxn.get_reaction_center()
   print(f"Reaction center: {reaction_center}")

Custom Similarity Metrics
~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

   from rxn_insight.reaction import Reaction
   from rxn_insight.utils import get_fp, get_similarity
   import numpy as np

   # Define two reactions
   rxn1 = "OB(O)c1ccccc1.Brc1ccccc1>>c1ccc(-c2ccccc2)cc1"
   rxn2 = "OB(O)c1ccc(C)cc1.Brc1ccccc1>>c1ccc(-c2ccc(C)cc2)cc1"

   # Get fingerprints
   fp1 = get_fp(rxn1, fp="Morgan", concatenate=True)
   fp2 = get_fp(rxn2, fp="Morgan", concatenate=True)

   # Calculate similarity using different metrics
   similarity_metrics = ["jaccard", "dice", "cosine", "euclidean", "manhattan"]

   for metric in similarity_metrics:
       similarity = get_similarity(fp1, fp2, metric=metric)
       print(f"{metric} similarity: {similarity:.4f}")

Working with Reaction Templates
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

   from rxn_insight.reaction import Reaction
   from rxn_insight.utils import get_reaction_template
   from rdkit import Chem
   from rdkit.Chem import AllChem

   # Create a reaction
   rxn_smiles = "OB(O)c1ccccc1.Brc1ccccc1>>c1ccc(-c2ccccc2)cc1"
   rxn = Reaction(rxn_smiles)

   # Extract template with different radii parameters
   template1 = get_reaction_template(rxn.mapped_reaction, radius_reactants=1, radius_products=1)
   template2 = get_reaction_template(rxn.mapped_reaction, radius_reactants=2, radius_products=1)

   print(f"Template (radius 1,1): {template1}")
   print(f"Template (radius 2,1): {template2}")

   # Use template to predict products for new reactants
   rxn_template = AllChem.ReactionFromSmarts(template1)
   new_reactants = ["OB(O)c1ccc(F)cc1", "Brc1ccc(Cl)cc1"]
   reactant_mols = [Chem.MolFromSmiles(r) for r in new_reactants]

   # Run the reaction
   products = rxn_template.RunReactants(reactant_mols)
   if products:
       predicted_product = Chem.MolToSmiles(products[0][0])
       print(f"Predicted product: {predicted_product}")

These examples demonstrate some of the advanced features available in
Rxn-INSIGHT. Refer to the source code for more detailed documentation of
each function and class.
