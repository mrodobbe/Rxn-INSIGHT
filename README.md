# Rxn-INSIGHT: Fast Chemical Reaction Analysis Using Bond-Electron Matrices

![Coverage Status](coverage-badge.svg)

Rxn-INSIGHT is an open-source algorithm, written in python, to classify and name chemical reactions, and suggest reaction conditions based on similarity and popularity.
* https://doi.org/10.1186/s13321-024-00834-z: Peer-reviewed publication on Rxn-INSIGHT
## 1. Installation
Rxn-INSIGHT relies on NumPy, Pandas, RDKit, RDChiral, and RXNMapper.

A virtual environment can be installed with Anaconda as follows:

```console
conda create -n rxn-insight python=3.10
conda activate rxn-insight
```

```
git clone https://github.com/schwallergroup/Rxn-INSIGHT.git
cd Rxn-INSIGHT
pip install .
```

Or, for developing with the optional dependencies, which are required to run the tests
and build the docs:
``` 
pip install -e ".[test,doc]"
```

All of the test environments can be run using the command `tox` from the top directory.
Alternatively, individual test environments can be run using the `-e` flag as 
in `tox -e env-name`. To run the tests, tests with coverage report, style checks, and
docs build, respectively:
```
tox -e py3
tox -e py3-coverage
tox -e style
tox -e docs
```

## 2. Usage

### Basic Usage
```python
from rxn_insight.reaction import Reaction
r = "c1ccccc1I.C=CC(=O)OC>>COC(=O)/C=C/c1ccccc1"  # Define a Reaction SMILES identifier
rxn = Reaction(r)
ri = rxn.get_reaction_info()
```

The reaction info contains most of the information:
```python
{'REACTION': 'C=CC(=O)OC.Ic1ccccc1>>COC(=O)/C=C/c1ccccc1', 
 'MAPPED_REACTION': '[CH3:1][O:2][C:3](=[O:4])[CH:5]=[CH2:6].I[c:7]1[cH:8][cH:9][cH:10][cH:11][cH:12]1>>[CH3:1][O:2][C:3](=[O:4])/[CH:5]=[CH:6]/[c:7]1[cH:8][cH:9][cH:10][cH:11][cH:12]1', 
 'N_REACTANTS': 2, 
 'N_PRODUCTS': 1, 
 'FG_REACTANTS': ('Aromatic halide', 'Vinyl'), 
 'FG_PRODUCTS': (), 
 'PARTICIPATING_RINGS_REACTANTS': ('c1ccccc1',), 
 'PARTICIPATING_RINGS_PRODUCTS': ('c1ccccc1',), 
 'ALL_RINGS_PRODUCTS': ('c1ccccc1',), 
 'BY-PRODUCTS': ('HI',), 
 'CLASS': 'C-C Coupling', 
 'TAG': '55becfded1a3842d5a03bbf3e1610411c659aff0806930400c4db2ef61f9c87f', 
 'SOLVENT': ('',), 
 'REAGENT': ('',), 
 'CATALYST': ('',), 
 'REF': '', 
 'NAME': 'Heck terminal vinyl', 
 'SCAFFOLD': 'c1ccccc1'}
```

### Similarity Search
A similarity search can be performed when a database with similar reactions is provided as a pandas DataFrame (df in this case). Another Pandas DataFrame is returned.
```python
df_nbs = rxn.find_neighbors(df, fp="MACCS", concatenate=True, threshold=0.5, broaden=True, full_search=False)
```

### Condition Suggestion
Reaction conditions can be suggested when a Pandas DataFrame is provided.
```python
rxn.suggest_conditions(df)
suggested_solvents = rxn.suggested_solvent
suggested_catalysts = rxn.suggested_catalyst
suggested_reagents = rxn.suggested_reagent
```

## 3. Datasets
The complete USPTO dataset that is analyzed by Rxn-INSIGHT, 
as described in the manuscript, can be found on 
Zenodo: https://doi.org/10.5281/zenodo.10171745. 
The `gzip` file should be downloaded and placed in the folder `data/`.

## 4. Reference
When using Rxn-INSIGHT for your own work, please refer to the original publication: <br>
`M. R. Dobbelaere, I. Lengyel, C. V. Stevens, and K. M. Van Geem, 
‘Rxn-INSIGHT: fast chemical reaction analysis using bond-electron matrices’, J. Cheminform., vol. 16, no. 1, Mar. 2024.`

```
@ARTICLE{Dobbelaere2024-es,
  title     = "{Rxn-INSIGHT}: fast chemical reaction analysis using
               bond-electron matrices",
  author    = "Dobbelaere, Maarten R and Lengyel, Istv{\'a}n and Stevens,
               Christian V and Van Geem, Kevin M",
  journal   = "J. Cheminform.",
  publisher = "Springer Science and Business Media LLC",
  volume    =  16,
  number    =  1,
  month     =  mar,
  year      =  2024,
  copyright = "https://creativecommons.org/licenses/by/4.0",
  language  = "en"
}
```
