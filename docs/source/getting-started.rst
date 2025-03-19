Getting Started
===============

This quick-start guide will get you up and running with Rxn-INSIGHT in
just a few minutes. We will cover installation, basic usage, and a simple
workflow for analyzing chemical reactions.

Step 1: Installation
--------------------

Install Rxn-INSIGHT using pip:

.. code:: bash

   pip install rxn_insight

Verify the installation by importing the package in Python:

.. code:: python

   import rxn_insight
   print("Rxn-INSIGHT installed successfully!")

Step 2: Your First Reaction Analysis
------------------------------------

Analyzing a simple Suzuki coupling reaction:

.. code:: python

   import rxn_insight as ri

   # Define a Suzuki coupling reaction in SMILES format
   reaction_smiles = "OB(O)c1ccccc1.Brc1ccccc1>>c1ccc(-c2ccccc2)cc1"

   # Create a Reaction object
   rxn = ri.Reaction(reaction_smiles)

   # Get reaction information
   info = rxn.get_reaction_info()

   # Print the basic information
   print(f"Reaction class: {info['CLASS']}")
   print(f"Reaction name: {info['NAME']}")
   print(f"Functional groups in reactants: {info['FG_REACTANTS']}")

You should see output similar to:

::

   Reaction class: C-C Coupling
   Reaction name: Suzuki coupling with boronic acids
   Functional groups in reactants: ['Aromatic halide', 'Boronic acid']

Step 3: Finding Similar Reactions
---------------------------------

To find similar reactions, you will need a reaction database. The USPTO
dataset processed for Rxn-INSIGHT can be downloaded from
`Zenodo <https://zenodo.org/records/10171745>`__.

.. code:: python

   import pandas as pd
   import rxn_insight as ri

   # Load the USPTO dataset
   # Download it first from https://zenodo.org/records/10171745
   df_uspto = pd.read_parquet("uspto_rxn_insight.gzip")

   # Define our reaction
   reaction_smiles = "OB(O)c1ccccc1.Brc1ccccc1>>c1ccc(-c2ccccc2)cc1"
   rxn = ri.Reaction(reaction_smiles)

   # Find similar reactions with at least 0.4 similarity
   similar_reactions = rxn.find_neighbors(df_uspto, threshold=0.4)

   # Print the top 3 most similar reactions
   for i, (idx, row) in enumerate(similar_reactions.head(3).iterrows()):
       print(f"\nSimilar reaction #{i+1} (Similarity: {row['SIMILARITY']:.2f}):")
       print(f"Reaction: {row['REACTION']}")
       print(f"Conditions: {row['SOLVENT']}, {row['CATALYST']}, {row['REAGENT']}")
       print(f"Yield: {row['YIELD']}%")

Step 4: Suggestion Reaction Conditions from Literature
------------------------------------------------------

Rxn-INSIGHT can recommend reaction conditions based on similar
reactions:

.. code:: python

   # Using the same reaction and database from Step 3
   conditions = rxn.suggest_conditions(df_uspto)

   print("\nRecommended conditions:")
   print(f"Solvent: {conditions['Solvent']}")
   print(f"Catalyst: {conditions['Catalyst']}")
   print(f"Reagent: {conditions['Reagent']}")

   # Get more detailed rankings
   print("\nTop 3 solvents:")
   solvent_ranking = rxn.suggested_solvent
   print(solvent_ranking[["NAME", "COUNT"]].head(3))

   print("\nTop 3 catalysts:")
   catalyst_ranking = rxn.suggested_catalyst
   print(catalyst_ranking[["NAME", "COUNT"]].head(3))

Bonus: Creating Your Own Reaction Database
------------------------------------------

If you have your own reaction data, you can create a custom database:

.. code:: python

   import rxn_insight as ri
   import pandas as pd

   # Create a simple DataFrame with your reactions
   data = {
       "reaction": [
           "OB(O)c1ccccc1.Brc1ccccc1>>c1ccc(-c2ccccc2)cc1",
           "CC(=O)c1ccccc1>>CC(O)c1ccccc1"
       ],
       "solvent": ["THF", "MeOH"],
       "reagent": ["K2CO3", "NaBH4"],
       "catalyst": ["Pd(PPh3)4", ""],
       "yield": [85, 92],
       "reference": ["Lab Notebook 1", "Lab Notebook 2"]
   }

   df = pd.DataFrame(data)

   # Create the database
   db = ri.Database()
   rxn_db = db.create_database_from_df(
       df,
       reaction_column="reaction",
       solvent_column="solvent",
       reagent_column="reagent",
       catalyst_column="catalyst",
       yield_column="yield",
       ref_column="reference"
   )

   # Save for future use
   db.save_to_parquet("my_reactions_database")

Next Steps
----------

Now that you have got the basics, you can:

1. Try analyzing different reactions
2. Build your own reaction database
3. Integrate condition prediction into your synthesis planning
4. Check out the detailed tutorials for advanced features


Happy reaction analyzing!
