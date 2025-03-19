import rxn_insight as ri


def test_load():
    rxn = ri.Reaction("CCCO>>CCC=O")


def test_molecule():
    m = ri.Molecule("CCCO")
