import rxn_insight as ri


def test_load():
    rxn = ri.Reaction("CCCO>>CCC=O")


def test_molecule():
    m = ri.Molecule("c1ccccc1", allow_pubchem=True)
    assert m.iupac_name.lower() == "benzene"


def test_molecule_rings():
    m = ri.Molecule("Cc1ccccc1", allow_pubchem=False)
    rings = m.get_rings()
    assert rings == ["c1ccccc1"]
