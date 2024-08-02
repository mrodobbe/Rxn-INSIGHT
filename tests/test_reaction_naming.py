from rxn_insight.reaction import Reaction


def test_protection():
    """
    Checks whether the protection reaction is named correctly.
    """
    rxn_smiles = "N[C@H](CO)Cc1ccccc1.CC(C)(C)OC(=O)OC(=O)O>>CC(C)(C)OC(=O)N[C@H](CO)Cc1ccccc1"
    rxn = Reaction(rxn_smiles)
    rxn_info = rxn.get_reaction_info()

    assert rxn_info["NAME"] == 'Boc amine protection of primary amine'


def test_acylation():
    """
    Checks whether the protection reaction is named correctly.
    """
    rxn_smiles = "O=C(Cl)CCl.CN>>CNC(=O)CCl"
    rxn = Reaction(rxn_smiles)
    rxn_info = rxn.get_reaction_info()

    assert rxn_info["CLASS"] == 'Acylation'
