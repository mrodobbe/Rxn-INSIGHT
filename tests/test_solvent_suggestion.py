import pandas as pd
from rxn_insight.reaction import Reaction
import os


def test_solvent_suggestion():
    test_dir = os.path.dirname(__file__)
    data_path = os.path.join(test_dir, '..', 'data', "example.gzip")
    df = pd.read_parquet(os.path.abspath(data_path))
    rxn_smiles = "BrC[P+](c1ccccc1)(c1ccccc1)c1ccccc1.CC1(C)OC2[C@H](O1)[C@@H](CC=O)O[C@H]2n1cnc2c(NC3CCCC3)ncnc21>>" \
                 "CC1(C)OC2[C@H](O1)[C@@H](C/C=C/Br)O[C@H]2n1cnc2c(NC3CCCC3)ncnc21"
    rxn = Reaction(rxn_smiles)
    rxn_info = rxn.get_reaction_info()
    rxn.suggest_conditions(df)
    rxn.find_neighbors(df, fp="morgan")
