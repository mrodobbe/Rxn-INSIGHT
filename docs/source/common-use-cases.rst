Potential Use Cases
=============================

This guide demonstrates how Rxn-INSIGHT can be applied to solve specific
chemistry problems.

Use Case 1: Reaction Classification and Analysis
------------------------------------------------

When planning a synthesis, it’s often useful to understand what type of
reaction you’re dealing with and its key components.

.. code:: python

   from rxn_insight.reaction import Reaction

   # Example 1: Analyzing a reduction reaction
   reduction_rxn = Reaction("CC(=O)c1ccccc1>>CC(O)c1ccccc1")
   info = reduction_rxn.get_reaction_info()

   print(f"Reaction class: {info['CLASS']}")
   print(f"Functional groups in reactants: {info['FG_REACTANTS']}")
   print(f"Functional groups in products: {info['FG_PRODUCTS']}")

Use Case 2: Reaction Condition Optimization
-------------------------------------------

When developing a new reaction or attempting to optimize an existing
one:

.. code:: python

   import pandas as pd
   from rxn_insight.reaction import Reaction

   # Load a reaction database
   df_rxns = pd.read_parquet("your_reaction_database.gzip")

   # Analyze a Buchwald-Hartwig amination reaction
   rxn = Reaction("Brc1ccccc1.NC1CCCCC1>>c1ccccc1NC1CCCCC1")

   # Get conditions from similar reactions
   conditions = rxn.suggest_conditions(df_rxns)

   # Get more detailed solvent rankings (sorted by frequency in similar reactions)
   solvent_ranking = rxn.suggested_solvent
   print(solvent_ranking.head(5))  # Top 5 solvents

   # Get more detailed catalyst rankings
   catalyst_ranking = rxn.suggested_catalyst
   print(catalyst_ranking.head(5))  # Top 5 catalysts

Use Case 3: Building a Reaction Database
----------------------------------------

Creating your own searchable reaction database:

.. code:: python

   import pandas as pd
   from rxn_insight.database import Database

   # Create a DataFrame with your reactions
   data = {
       "reaction": [
           "OB(O)c1ccccc1.Brc1ccccc1>>c1ccc(-c2ccccc2)cc1",
           "CC(=O)c1ccccc1>>CC(O)c1ccccc1",
           "Brc1ccccc1.NC1CCCCC1>>c1ccccc1NC1CCCCC1"
       ],
       "solvent": ["THF", "MeOH", "Toluene"], 
       "catalyst": ["Pd(PPh3)4", "", "Pd2(dba)3"],
       "reagent": ["K2CO3", "NaBH4", "t-BuONa"],
       "yield": [85, 92, 78],
       "reference": ["Paper1", "Paper2", "Paper3"]
   }

   df = pd.DataFrame(data)

   # Initialize Database object
   db = Database()

   # Create the database with analyzed reactions
   rxn_db = db.create_database_from_df(
       df,
       reaction_column="reaction",
       solvent_column="solvent",
       reagent_column="reagent",
       catalyst_column="catalyst",
       yield_column="yield",
       ref_column="reference"
   )

   # Save the database for future use
   db.save_to_parquet("my_reaction_database")

   # View class distribution in your database
   class_distribution = db.get_class_distribution()
   print(class_distribution)

Use Case 4: Compound Analysis
-----------------------------

Identifying and working with molecular scaffolds:

.. code:: python

   from rxn_insight.reaction import Molecule

   # Create a Molecule object from SMILES
   mol = Molecule("c1ccc(-c2ccccc2)cc1")  # Biphenyl

   # Get the scaffold
   scaffold = mol.scaffold
   print(f"Scaffold: {scaffold}")

   # Get functional groups
   functional_groups = mol.get_functional_groups()
   print(f"Functional groups: {functional_groups}")

   # Get ring systems
   rings = mol.get_rings()
   print(f"Rings: {rings}")

   # Find reactions that produce this product
   # (if you have a reaction database)
   import pandas as pd
   df_rxns = pd.read_parquet("your_reaction_database.gzip")
   related_reactions = mol.search_reactions(df_rxns)

Use Case 5: Reaction Template Extraction
----------------------------------------

Extract generic reaction templates that can be applied to similar
compounds:

.. code:: python

   from rxn_insight.reaction import Reaction
   from rxn_insight.utils import get_reaction_template

   # Define a specific reaction
   reaction_smiles = "OB(O)c1ccccc1.Brc1ccccc1>>c1ccc(-c2ccccc2)cc1"
   rxn = Reaction(reaction_smiles)

   # Get the reaction template
   # The template will contain the core reaction transformation
   # with generalized atom environments
   template = rxn.template

   print(f"Reaction template: {template}")

   # You can also extract a template directly from a mapped reaction
   from rxnmapper import RXNMapper
   rxn_mapper = RXNMapper()
   mapped_rxn = rxn_mapper.get_attention_guided_atom_maps([reaction_smiles])[0]["mapped_rxn"]

   # Extract template with custom radius parameters
   template = get_reaction_template(mapped_rxn, radius_reactants=2, radius_products=1)
