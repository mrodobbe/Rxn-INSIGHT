# Rxn-INSIGHT: Fast Chemical Reaction Analysis Using Bond-Electron Matrices

![Coverage Status](coverage-badge.svg)

Rxn-INSIGHT is an open-source algorithm, written in python, to classify and name chemical reactions, and suggest reaction conditions based on similarity and popularity.
* https://doi.org/10.1186/s13321-024-00834-z: Peer-reviewed publication on Rxn-INSIGHT

> **See also: [ReactionClassifier](https://github.com/schwallergroup/ReactionClassifier)** — our
> newer reaction naming tool, which uses agentically generated, verifiable rules for
> deterministic and self-expanding reaction classification.
>
> D. Armstrong, M. R. Dobbelaere, V. Olikauskas, H. Avila, O. Susanu, J. Waser, and P. Schwaller,
> *Agentic generation of verifiable rules for deterministic, self-expanding reaction
> classification*, [arXiv:2607.01061](https://arxiv.org/abs/2607.01061).

## 1. Installation
Rxn-INSIGHT relies on NumPy, SciPy, Pandas, RDKit, RDChiral, joblib, and RXNMapper.
It requires Python 3.10 or newer and is tested on 3.10 through 3.14.

A virtual environment can be installed with Anaconda as follows:

```console
conda create -n rxn-insight python=3.11
conda activate rxn-insight
```

### Option 1: Installing via PyPI:
```
pip install rxn-insight
```

Support for the [Open Reaction Database](https://open-reaction-database.org) is optional:
```
pip install "rxn-insight[ord]"
```
### Option 2: Installing directly from source:
```
git clone https://github.com/mrodobbe/Rxn-INSIGHT.git
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
 'FG_REACTANTS': ['Aromatic halide', 'Vinyl'], 
 'FG_PRODUCTS': [], 
 'PARTICIPATING_RINGS_REACTANTS': ['c1ccccc1'], 
 'PARTICIPATING_RINGS_PRODUCTS': ['c1ccccc1'], 
 'ALL_RINGS_PRODUCTS': ['c1ccccc1'], 
 'BY-PRODUCTS': ['HI'], 
 'CLASS': 'C-C Coupling', 
 'TAG': '55becfded1a3842d5a03bbf3e1610411c659aff0806930400c4db2ef61f9c87f', 
 'SOLVENT': [''], 
 'REAGENT': [''], 
 'CATALYST': [''], 
 'REF': '', 
 'NAME': 'Heck terminal vinyl', 
 'SCAFFOLD': 'c1ccccc1'}
```

### Explaining a Reaction
`explain()` reports *what changed* rather than just how the reaction is labelled —
the bonds formed and broken, the atoms lost, the rings that survived, and the
template that describes the transformation:

```python
rxn = Reaction("c1ccccc1I.C=CC(=O)OC>>COC(=O)/C=C/c1ccccc1")
explanation = rxn.explain()

for change in explanation["bond_changes"]:
    print(change)
# C-C single bond formed between atoms 6 and 7.
# C-I single bond broken between atoms 7 and 13.

explanation["leaving_groups"]   # ['Atom types I, with atom indices 12']
explanation["rings"]            # {'preserved': ['benzene']}
explanation["classification"]   # {'class': 'C-C Coupling', 'name': 'Heck terminal vinyl'}
```

The full set of keys is `reaction_center`, `bond_changes`,
`stereochemistry_changes`, `leaving_groups`, `added_atoms`, `rings`,
`functional_groups`, `template`, and `classification`.

### Reaction Naming
If you only need the *name* of a reaction, the top-level `name_reaction` helper is the fastest route.
It matches a reaction SMILES against the SMIRKS database directly — without building a full
`Reaction` object or atom-mapping — so it is far faster than the full pipeline (naming is fully
decoupled from classification):

```python
import rxn_insight as ri

ri.name_reaction("c1ccccc1I.C=CC(=O)OC>>COC(=O)/C=C/c1ccccc1")
# 'Heck terminal vinyl'
```

To name many reactions at once, `batch_name_reaction` runs in parallel across `n_jobs` worker
processes and returns a list of names (`"OtherReaction"` when nothing matches):

```python
names = ri.batch_name_reaction(
    ["c1ccccc1I.C=CC(=O)OC>>COC(=O)/C=C/c1ccccc1", "CCO>>CC=O"],
    n_jobs=8,
)
# ['Heck terminal vinyl', 'Oxidation or Dehydrogenation of Alcohols to Aldehydes and Ketones']
```

For a rule-based classifier that expands its own rule set, see
[ReactionClassifier](https://github.com/schwallergroup/ReactionClassifier)
([arXiv:2607.01061](https://arxiv.org/abs/2607.01061)).

### Looking Up Compounds by Name
`Compound` resolves a chemical name to a structure (via OPSIN by default, or PubChem)
and then behaves like a `Molecule`:

```python
from rxn_insight.molecule import Compound

benzene = Compound("benzene")
benzene.smiles      # 'c1ccccc1'
benzene.scaffold    # 'c1ccccc1'
```

Pass `use_opsin=False` to resolve through PubChem instead, or `allow_pubchem=True`
to additionally fetch descriptive information.

Both routes need network access. A name that cannot be resolved raises `ValueError`
on either route, so a `Compound` that is constructed successfully is always fully
initialised.

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

### Creating a Rxn-INSIGHT-Compatible Database
If you want to use similarity search or condition suggestion with your own data, then you need to create a database 
that is compatible with Rxn-INSIGHT. This action is possible with the Database module from either a csv file 
or a Pandas DataFrame. Below an example is shown:

```python
from rxn_insight.database import Database
from rxn_insight.reaction import Reaction

db = Database()
db.create_database_from_csv(fname="test.csv", 
                            reaction_column="RXN", 
                            solvent_column="SOLVENTS")
r = Reaction("CCO>>CC=O")
sug_conds = r.suggest_conditions(db.df)
```

Two arguments are required: `fname` (location of your csv file) and 
`reaction_column` (the name of the column containing the reactions).
The optional arguments are `solvent_column`, `reagent_column`, `catalyst_column`, `yield_column`, and `ref_column`. 
Default (or when these columns are not available), the column names will be set to respectively SOLVENT, REAGENT, CATALYST, YIELD, REF.

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
