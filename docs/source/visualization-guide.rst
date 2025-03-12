Visualizing Chemical Reactions
==============================

This tutorial demonstrates how to create effective visualizations of
reactions, their components, and analysis results using Rxn-INSIGHT
combined with common Python visualization libraries.

Prerequisites
-------------

- Rxn-INSIGHT installed
- Basic knowledge of Python
- Required packages:

  - matplotlib
  - seaborn
  - RDKit
  - pandas
  - numpy

Install the visualization packages:

.. code:: bash

   pip install matplotlib seaborn

1. Visualizing Chemical Structures
----------------------------------

First, let’s see how to visualize the molecules in a reaction:

.. code:: python

   from rxn_insight.reaction import Reaction
   from rdkit import Chem
   from rdkit.Chem import Draw
   import matplotlib.pyplot as plt

   # Define a reaction
   reaction_smiles = "OB(O)c1ccccc1.Brc1ccccc1>>c1ccc(-c2ccccc2)cc1"
   rxn = Reaction(reaction_smiles)

   # Extract reactants and product
   reactants_smiles = rxn.reactants.split(".")
   product_smiles = rxn.products

   # Convert to RDKit molecules
   reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]
   product_mol = Chem.MolFromSmiles(product_smiles)

   # Create a visualization of reactants and products
   fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))

   # Visualize reactants
   reactant_img = Draw.MolsToGridImage(reactant_mols, molsPerRow=2, subImgSize=(200, 200))
   ax1.imshow(reactant_img)
   ax1.set_title("Reactants")
   ax1.axis("off")

   # Visualize product
   product_img = Draw.MolToImage(product_mol, size=(200, 200))
   ax2.imshow(product_img)
   ax2.set_title("Product")
   ax2.axis("off")

   plt.tight_layout()
   plt.savefig('reaction_molecules.png', dpi=300)
   plt.show()

2. Visualizing Reaction Center
------------------------------

Highlight the reaction center to focus on the transforming parts:

.. code:: python

   from rxn_insight.reaction import Reaction
   from rxn_insight.utils import draw_chemical_reaction

   # Define a reaction
   reaction_smiles = "OB(O)c1ccccc1.Brc1ccccc1>>c1ccc(-c2ccccc2)cc1"
   rxn = Reaction(reaction_smiles)

   # Get the reaction info to identify the mapped reaction
   reaction_info = rxn.get_reaction_info()
   mapped_rxn = reaction_info["MAPPED_REACTION"]

   # Draw the reaction with highlighting
   reaction_svg = draw_chemical_reaction(mapped_rxn, highlightByReactant=True)

   # Save the SVG to a file
   with open("reaction_center.svg", "w") as f:
       f.write(reaction_svg)

   # To display in a Jupyter notebook:
   from IPython.display import SVG
   SVG(reaction_svg)

3. Visualizing Reaction Classification Results
----------------------------------------------

Create a pie chart of reaction classes in your database:

.. code:: python

   import pandas as pd
   import matplotlib.pyplot as plt
   import seaborn as sns
   from rxn_insight.database import Database

   # Load or create a database
   db = Database()
   try:
       # Try to load an existing database
       df = pd.read_parquet("my_rxn_database.gzip")
       db.df = df
   except:
       # Or use the USPTO dataset if available
       df = pd.read_parquet("uspto_rxn_insight.gzip")
       db.df = df

   # Get the class distribution
   class_distribution = db.get_class_distribution()

   # Create a pie chart
   plt.figure(figsize=(10, 8))
   plt.pie(class_distribution["COUNT"], 
           labels=class_distribution["CLASS"], 
           autopct='%1.1f%%',
           textprops={'fontsize': 9})
   plt.title("Distribution of Reaction Classes", fontsize=14)
   plt.tight_layout()
   plt.savefig('reaction_classes.png', dpi=300)
   plt.show()

4. Visualizing Similarity Search Results
----------------------------------------

Create a heatmap of similarity between various reactions:

.. code:: python

   import numpy as np
   import pandas as pd
   import matplotlib.pyplot as plt
   import seaborn as sns
   from rxn_insight.reaction import Reaction
   from rxn_insight.utils import get_fp, get_similarity

   # Define several reactions
   reaction_smiles = [
       "OB(O)c1ccccc1.Brc1ccccc1>>c1ccc(-c2ccccc2)cc1",  # Suzuki coupling
       "OB(O)c1ccc(F)cc1.Brc1ccc(Cl)cc1>>c1ccc(F)c(-c2ccc(Cl)cc2)c1",  # Similar Suzuki
       "CC(=O)c1ccccc1>>CC(O)c1ccccc1",  # Reduction
       "CC(O)c1ccccc1>>CC(=O)c1ccccc1",  # Oxidation
       "Brc1ccccc1.NC1CCCCC1>>c1ccccc1NC1CCCCC1"  # Buchwald-Hartwig
   ]

   # Get fingerprints
   fps = [get_fp(rxn, fp="MACCS", concatenate=True) for rxn in reaction_smiles]

   # Calculate pairwise similarities
   n = len(fps)
   similarity_matrix = np.zeros((n, n))

   for i in range(n):
       for j in range(n):
           similarity_matrix[i, j] = get_similarity(fps[i], fps[j], metric="jaccard")

   # Create labels for the reactions
   labels = [
       "Suzuki Coupling (basic)",
       "Suzuki Coupling (substituted)",
       "Ketone Reduction",
       "Alcohol Oxidation",
       "Buchwald-Hartwig Amination"
   ]

   # Create a heatmap
   plt.figure(figsize=(10, 8))
   sns.heatmap(similarity_matrix, annot=True, fmt=".2f", cmap="YlGnBu",
               xticklabels=labels, yticklabels=labels)
   plt.title("Reaction Similarity Matrix (Jaccard Index)")
   plt.tight_layout()
   plt.savefig('reaction_similarity.png', dpi=300)
   plt.show()

5. Visualizing Condition Recommendations
----------------------------------------

Create bar charts of recommended solvents, catalysts, and reagents:

.. code:: python

   import pandas as pd
   import matplotlib.pyplot as plt
   import seaborn as sns
   from rxn_insight.reaction import Reaction

   # Define a reaction
   reaction_smiles = "OB(O)c1ccccc1.Brc1ccccc1>>c1ccc(-c2ccccc2)cc1"
   rxn = Reaction(reaction_smiles)

   # Load a reaction database
   df_rxns = pd.read_parquet("uspto_rxn_insight.gzip")

   # Get condition suggestions
   conditions = rxn.suggest_conditions(df_rxns)

   # Get detailed rankings
   solvent_ranking = rxn.suggested_solvent
   catalyst_ranking = rxn.suggested_catalyst
   reagent_ranking = rxn.suggested_reagent

   # Create a figure with 3 subplots
   fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 15))

   # Plot top 5 solvents
   top_solvents = solvent_ranking.head(5)
   sns.barplot(x='COUNT', y='NAME', data=top_solvents, ax=ax1, palette='Blues_d')
   ax1.set_title('Top 5 Recommended Solvents')
   ax1.set_xlabel('Count')
   ax1.set_ylabel('Solvent')

   # Plot top 5 catalysts
   top_catalysts = catalyst_ranking.head(5)
   sns.barplot(x='COUNT', y='NAME', data=top_catalysts, ax=ax2, palette='Greens_d')
   ax2.set_title('Top 5 Recommended Catalysts')
   ax2.set_xlabel('Count')
   ax2.set_ylabel('Catalyst')

   # Plot top 5 reagents
   top_reagents = reagent_ranking.head(5)
   sns.barplot(x='COUNT', y='NAME', data=top_reagents, ax=ax3, palette='Reds_d')
   ax3.set_title('Top 5 Recommended Reagents')
   ax3.set_xlabel('Count')
   ax3.set_ylabel('Reagent')

   plt.tight_layout()
   plt.savefig('reaction_conditions.png', dpi=300)
   plt.show()

6. Visualizing Reaction Networks
--------------------------------

Create a graph showing relationships between reactions:

.. code:: python

   import pandas as pd
   import matplotlib.pyplot as plt
   import networkx as nx
   from rxn_insight.reaction import Reaction, Molecule

   # Install networkx if you don't have it
   # pip install networkx

   # Define a target product
   target_smiles = "c1ccc(-c2ccccc2)cc1"  # Biphenyl
   target = Molecule(target_smiles)

   # Load a reaction database
   df_rxns = pd.read_parquet("uspto_rxn_insight.gzip")

   # Find reactions that produce this target
   target_reactions = target.search_reactions(df_rxns)

   if target_reactions is not None and len(target_reactions) > 0:
       # Create a graph
       G = nx.DiGraph()
       
       # Add the target node
       G.add_node(target_smiles, type="product", label="Target")
       
       # Add reactions and reactants
       for i, row in target_reactions.head(5).iterrows():
           # Add reaction node
           reaction_id = f"Reaction {i}"
           G.add_node(reaction_id, type="reaction", 
                      label=f"Class: {row['CLASS']}", yield=row['YIELD'])
           
           # Add edge from reaction to product
           G.add_edge(reaction_id, target_smiles)
           
           # Add reactants
           reactants = row['REACTION'].split(">>")[0].split(".")
           for j, reactant in enumerate(reactants):
               reactant_id = f"{reactant}_{i}_{j}"
               G.add_node(reactant_id, type="reactant", label=f"Reactant {j+1}")
               G.add_edge(reactant_id, reaction_id)
       
       # Create positions for the graph
       pos = nx.spring_layout(G, k=0.5, iterations=50)
       
       # Draw the graph
       plt.figure(figsize=(12, 10))
       
       # Draw nodes by type
       product_nodes = [n for n, d in G.nodes(data=True) if d.get('type') == 'product']
       reaction_nodes = [n for n, d in G.nodes(data=True) if d.get('type') == 'reaction']
       reactant_nodes = [n for n, d in G.nodes(data=True) if d.get('type') == 'reactant']
       
       nx.draw_networkx_nodes(G, pos, nodelist=product_nodes, node_color='red', 
                             node_size=500, alpha=0.8)
       nx.draw_networkx_nodes(G, pos, nodelist=reaction_nodes, node_color='blue', 
                             node_size=400, alpha=0.8)
       nx.draw_networkx_nodes(G, pos, nodelist=reactant_nodes, node_color='green', 
                             node_size=300, alpha=0.8)
       
       # Draw edges
       nx.draw_networkx_edges(G, pos, width=1.0, alpha=0.5, arrows=True)
       
       # Draw labels
       labels = {n: G.nodes[n].get('label', n) for n in G.nodes()}
       # Simplify SMILES labels
       for n in labels:
           if isinstance(n, str) and '(' in n and len(n) > 15:
               labels[n] = "Molecule"
       
       nx.draw_networkx_labels(G, pos, labels=labels, font_size=8)
       
       plt.title("Reaction Network for Target Product", fontsize=14)
       plt.axis('off')
       plt.tight_layout()
       plt.savefig('reaction_network.png', dpi=300)
       plt.show()

7. Interactive Visualization with Plotly
----------------------------------------

For interactive visualizations (especially useful in Jupyter notebooks):

.. code:: python

   import pandas as pd
   import plotly.express as px
   import plotly.graph_objects as go
   from rxn_insight.database import Database

   # Install plotly if you don't have it
   # pip install plotly

   # Load or create a database
   db = Database()
   df = pd.read_parquet("uspto_rxn_insight.gzip")
   db.df = df

   # Get class distribution
   class_dist = db.get_class_distribution()

   # Create an interactive pie chart
   fig = px.pie(class_dist, values='COUNT', names='CLASS', 
                title='Interactive Distribution of Reaction Classes')
   fig.update_traces(textposition='inside', textinfo='percent+label')
   fig.show()

   # Interactive bar chart for reaction yields by class
   # First, get average yields by class
   yield_by_class = df.groupby('CLASS')['YIELD'].mean().reset_index()

   # Create interactive bar chart
   fig = px.bar(yield_by_class, x='CLASS', y='YIELD',
                title='Average Yield by Reaction Class',
                labels={'YIELD': 'Average Yield (%)', 'CLASS': 'Reaction Class'},
                color='YIELD', color_continuous_scale=px.colors.sequential.Viridis)
   fig.update_layout(xaxis_tickangle=-45)
   fig.show()

   # Save interactive visualization as HTML
   fig.write_html("interactive_reaction_stats.html")

8. Custom Visualization Functions
---------------------------------

Create reusable functions for common visualizations:

.. code:: python

   import pandas as pd
   import matplotlib.pyplot as plt
   import seaborn as sns
   from rdkit import Chem
   from rdkit.Chem import Draw
   from rxn_insight.reaction import Reaction

   def visualize_reaction(reaction_smiles, title=None, save_path=None):
       """Create a visualization of a reaction with reactants, arrow, and products."""
       rxn = Reaction(reaction_smiles)
       
       # Extract reactants and product
       reactants_smiles = rxn.reactants.split(".")
       product_smiles = rxn.products.split(".")
       
       # Convert to RDKit molecules
       reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]
       product_mols = [Chem.MolFromSmiles(p) for p in product_smiles]
       
       # Create the reaction drawing
       rxn_drawing = Draw.ReactionToImage(reactant_mols, product_mols)
       
       # Create figure
       plt.figure(figsize=(10, 4))
       plt.imshow(rxn_drawing)
       
       if title:
           plt.title(title)
       plt.axis('off')
       
       if save_path:
           plt.savefig(save_path, dpi=300, bbox_inches='tight')
       
       plt.show()

   def plot_condition_distribution(df, condition_col, title=None, top_n=10, save_path=None):
       """Plot the distribution of a condition (solvent, catalyst, reagent) in a database."""
       # Get counts for each unique value
       condition_counts = df[condition_col].value_counts().reset_index()
       condition_counts.columns = [condition_col, 'COUNT']
       
       # Take the top N most common
       top_conditions = condition_counts.head(top_n)
       
       # Create plot
       plt.figure(figsize=(10, 6))
       sns.barplot(x='COUNT', y=condition_col, data=top_conditions)
       
       if title:
           plt.title(title)
       else:
           plt.title(f'Distribution of {condition_col}')
       
       plt.tight_layout()
       
       if save_path:
           plt.savefig(save_path, dpi=300)
       
       plt.show()
       
       return top_conditions

   # Example usage
   reaction_smiles = "OB(O)c1ccccc1.Brc1ccccc1>>c1ccc(-c2ccccc2)cc1"
   visualize_reaction(reaction_smiles, title="Suzuki Coupling", save_path="suzuki.png")

   # Load a database and plot condition distributions
   df = pd.read_parquet("uspto_rxn_insight.gzip")
   plot_condition_distribution(df, "SOLVENT", title="Most Common Solvents", save_path="solvents.png")

9. Visualizing Functional Groups and Rings
------------------------------------------

Create visualizations to highlight functional groups and ring systems:

.. code:: python

   from rxn_insight.reaction import Reaction, Molecule
   from rdkit import Chem
   from rdkit.Chem import Draw, AllChem
   import matplotlib.pyplot as plt

   def highlight_functional_groups(mol_smiles, title=None, save_path=None):
       """Highlight functional groups in a molecule."""
       # Create a Molecule object to get functional groups
       mol_obj = Molecule(mol_smiles)
       fg_list = mol_obj.get_functional_groups()
       
       # Create RDKit molecule
       mol = Chem.MolFromSmiles(mol_smiles)
       
       # Create a list to hold highlighted versions
       highlighted_mols = []
       legends = []
       
       # Create the original molecule with no highlighting
       highlighted_mols.append(mol)
       legends.append("Original")
       
       # Highlight each functional group
       for fg in fg_list:
           # Create a copy of the molecule
           mol_copy = Chem.Mol(mol)
           
           # Try to find the functional group with a SMARTS pattern
           try:
               # Common functional group SMARTS patterns
               fg_smarts = {
                   "Alcohol": "[OX2H]",
                   "Aldehyde": "[CX3H1](=O)[#6]",
                   "Ketone": "[#6][CX3](=O)[#6]",
                   "Carboxylic acid": "[CX3](=O)[OX2H1]",
                   "Ester": "[#6][CX3](=O)[OX2][#6]",
                   "Amide": "[NX3][CX3](=[OX1])[#6]",
                   "Amine": "[NX3;H2,H1,H0;!$(NC=O)]",
                   "Nitrile": "[NX1]#[CX2]",
                   "Nitro": "[$([NX3](=O)=O),$([NX3+](=O)[O-])]",
                   "Aromatic": "c1ccccc1",
                   "Boronic acid": "[B][O][H]",
                   "Aromatic halide": "[c][F,Cl,Br,I]",
                   # Add more patterns as needed
               }
               
               if fg in fg_smarts:
                   pattern = Chem.MolFromSmarts(fg_smarts[fg])
                   if pattern:
                       matches = mol_copy.GetSubstructMatches(pattern)
                       if matches:
                           # Combine all matches into one highlight
                           atoms_to_highlight = set()
                           for match in matches:
                               atoms_to_highlight.update(match)
                           
                           # Highlight the atoms
                           for atom_idx in atoms_to_highlight:
                               mol_copy.GetAtomWithIdx(atom_idx).SetProp("atomNote", fg)
                           
                           highlighted_mols.append(mol_copy)
                           legends.append(fg)
           except Exception as e:
               print(f"Error highlighting {fg}: {e}")
       
       # Create the grid image
       grid_img = Draw.MolsToGridImage(
           highlighted_mols, 
           molsPerRow=min(3, len(highlighted_mols)),
           subImgSize=(250, 250),
           legends=legends,
           useSVG=False
       )
       
       # Display the image
       plt.figure(figsize=(12, 10))
       plt.imshow(grid_img)
       
       if title:
           plt.title(title)
       plt.axis('off')
       
       if save_path:
           plt.savefig(save_path, dpi=300, bbox_inches='tight')
       
       plt.show()

   def highlight_rings(mol_smiles, title=None, save_path=None):
       """Highlight ring systems in a molecule."""
       # Create a Molecule object to get rings
       mol_obj = Molecule(mol_smiles)
       ring_list = mol_obj.get_rings()
       
       # Create RDKit molecule
       mol = Chem.MolFromSmiles(mol_smiles)
       
       # Create a list to hold highlighted versions
       highlighted_mols = []
       legends = []
       
       # Create the original molecule with no highlighting
       highlighted_mols.append(mol)
       legends.append("Original")
       
       # Highlight each ring system
       for i, ring_smiles in enumerate(ring_list):
           # Create a copy of the molecule
           mol_copy = Chem.Mol(mol)
           
           try:
               # Convert ring SMILES to a substructure
               ring_mol = Chem.MolFromSmiles(ring_smiles)
               if ring_mol:
                   # Generate 3D coordinates to help with matching
                   AllChem.EmbedMolecule(ring_mol)
                   AllChem.EmbedMolecule(mol_copy)
                   
                   # Find matches
                   matches = mol_copy.GetSubstructMatches(ring_mol)
                   if matches:
                       # Combine all matches
                       atoms_to_highlight = set()
                       for match in matches:
                           atoms_to_highlight.update(match)
                       
                       # Highlight the atoms
                       for atom_idx in atoms_to_highlight:
                           mol_copy.GetAtomWithIdx(atom_idx).SetProp("atomNote", f"Ring {i+1}")
                       
                       highlighted_mols.append(mol_copy)
                       legends.append(f"Ring {i+1}: {ring_smiles}")
           except Exception as e:
               print(f"Error highlighting ring {i}: {e}")
       
       # Create the grid image
       grid_img = Draw.MolsToGridImage(
           highlighted_mols, 
           molsPerRow=min(2, len(highlighted_mols)),
           subImgSize=(300, 300),
           legends=legends,
           useSVG=False
       )
       
       # Display the image
       plt.figure(figsize=(12, 10))
       plt.imshow(grid_img)
       
       if title:
           plt.title(title)
       plt.axis('off')
       
       if save_path:
           plt.savefig(save_path, dpi=300, bbox_inches='tight')
       
       plt.show()

   # Example usage
   mol_smiles = "c1ccc(-c2ccccc2)cc1"  # Biphenyl
   highlight_functional_groups(mol_smiles, title="Functional Groups in Biphenyl", save_path="biphenyl_fg.png")
   highlight_rings(mol_smiles, title="Ring Systems in Biphenyl", save_path="biphenyl_rings.png")

10. Dashboard for Reaction Analysis
-----------------------------------

Finally, let’s create a simple dashboard-style view that integrates
multiple visualizations:

.. code:: python

   import pandas as pd
   import matplotlib.pyplot as plt
   from matplotlib.gridspec import GridSpec
   import seaborn as sns
   from rxn_insight.reaction import Reaction
   from rdkit import Chem
   from rdkit.Chem import Draw

   def create_reaction_dashboard(reaction_smiles, database_path=None, save_path=None):
       """Create a comprehensive dashboard for a reaction."""
       # Initialize the reaction
       rxn = Reaction(reaction_smiles)
       info = rxn.get_reaction_info()
       
       # Set up the figure
       fig = plt.figure(figsize=(15, 12))
       gs = GridSpec(3, 3, figure=fig)
       
       # 1. Basic reaction information
       ax_info = fig.add_subplot(gs[0, 0:2])
       ax_info.axis('off')
       info_text = (
           f"Reaction Class: {info['CLASS']}\n"
           f"Reaction Name: {info['NAME']}\n"
           f"Reactants: {info['N_REACTANTS']}\n"
           f"Products: {info['N_PRODUCTS']}\n"
           f"Functional Groups in Reactants: {', '.join(info['FG_REACTANTS'])}\n"
           f"By-products: {', '.join(info['BY-PRODUCTS'])}\n"
           f"Scaffold: {info['SCAFFOLD']}"
       )
       ax_info.text(0.05, 0.95, info_text, transform=ax_info.transAxes, 
                   verticalalignment='top', fontsize=10)
       ax_info.set_title("Reaction Information", fontsize=12)
       
       # 2. Reaction visualization
       ax_rxn = fig.add_subplot(gs[0, 2])
       # Extract reactants and products
       reactants_smiles = rxn.reactants.split(".")
       products_smiles = rxn.products.split(".")
       # Convert to RDKit molecules
       reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]
       product_mols = [Chem.MolFromSmiles(p) for p in products_smiles]
       # Draw the reaction
       rxn_img = Draw.ReactionToImage(reactant_mols, product_mols, subImgSize=(200, 150))
       ax_rxn.imshow(rxn_img)
       ax_rxn.axis('off')
       ax_rxn.set_title("Reaction Visualization", fontsize=12)
       
       # Load database if available
       if database_path:
           try:
               df_rxns = pd.read_parquet(database_path)
               
               # 3. Condition recommendations
               conditions = rxn.suggest_conditions(df_rxns)
               
               # Solvents
               ax_solvent = fig.add_subplot(gs[1, 0])
               solvent_ranking = rxn.suggested_solvent.head(5)
               sns.barplot(x='COUNT', y='NAME', data=solvent_ranking, ax=ax_solvent)
               ax_solvent.set_title("Recommended Solvents", fontsize=12)
               
               # Catalysts
               ax_catalyst = fig.add_subplot(gs[1, 1])
               catalyst_ranking = rxn.suggested_catalyst.head(5)
               sns.barplot(x='COUNT', y='NAME', data=catalyst_ranking, ax=ax_catalyst)
               ax_catalyst.set_title("Recommended Catalysts", fontsize=12)
               
               # Reagents
               ax_reagent = fig.add_subplot(gs[1, 2])
               reagent_ranking = rxn.suggested_reagent.head(5)
               sns.barplot(x='COUNT', y='NAME', data=reagent_ranking, ax=ax_reagent)
               ax_reagent.set_title("Recommended Reagents", fontsize=12)
               
               # 4. Similar reactions
               ax_similar = fig.add_subplot(gs[2, 0:3])
               similar_reactions = rxn.find_neighbors(df_rxns, threshold=0.3, max_return=5)
               
               if similar_reactions is not None and len(similar_reactions) > 0:
                   similar_text = "Top Similar Reactions:\n\n"
                   for i, (idx, row) in enumerate(similar_reactions.head(5).iterrows()):
                       similar_text += (
                           f"{i+1}. Similarity: {row['SIMILARITY']:.2f}\n"
                           f"   Reaction: {row['REACTION']}\n"
                           f"   Conditions: Solvent: {row['SOLVENT']}, "
                           f"Catalyst: {row['CATALYST']}, Reagent: {row['REAGENT']}\n"
                           f"   Yield: {row['YIELD']}%\n\n"
                       )
               else:
                   similar_text = "No similar reactions found in the database."
                   
               ax_similar.axis('off')
               ax_similar.text(0.05, 0.95, similar_text, transform=ax_similar.transAxes, 
                             verticalalignment='top', fontsize=9)
               ax_similar.set_title("Similar Reactions", fontsize=12)
               
           except Exception as e:
               ax_error = fig.add_subplot(gs[1:, :])
               ax_error.axis('off')
               ax_error.text(0.5, 0.5, f"Error loading or processing database: {e}", 
                            transform=ax_error.transAxes, 
                            horizontalalignment='center', 
                            verticalalignment='center',
                            fontsize=12, color='red')
       else:
           ax_note = fig.add_subplot(gs[1:, :])
           ax_note.axis('off')
           ax_note.text(0.5, 0.5, "No database provided. Condition recommendations not available.", 
                      transform=ax_note.transAxes, 
                      horizontalalignment='center', 
                      verticalalignment='center',
                      fontsize=12)
       
       # Add title to the figure
       fig.suptitle(f"Comprehensive Analysis of {info['NAME']}", fontsize=16, y=0.98)
       
       plt.tight_layout()
       
       if save_path:
           plt.savefig(save_path, dpi=300, bbox_inches='tight')
       
       plt.show()

   # Example usage
   reaction_smiles = "OB(O)c1ccccc1.Brc1ccccc1>>c1ccc(-c2ccccc2)cc1"
   create_reaction_dashboard(reaction_smiles, 
                            database_path="uspto_rxn_insight.gzip", 
                            save_path="reaction_dashboard.png")

This dashboard provides a comprehensive overview of the reaction,
combining information about its classification, structure, recommended
conditions, and similar reactions from the database.

By using these visualization techniques, you can better communicate the
insights gained from Rxn-INSIGHT’s analysis and make more informed
decisions for reaction optimization and synthesis planning.
