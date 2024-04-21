from rxn_insight.classification import ReactionClassifier


def test_initialization():
    rxn_smiles_with_atom_mapping = "[CH3:9][CH:8]([CH3:10])[c:7]1[cH:11][cH:12][cH:13][cH:14][cH:15]1.[CH2:3]([CH2:4][C:5](=[O:6])Cl)[CH2:2][Cl:1].[Al+3].[Cl-].[Cl-].[Cl-].C(Cl)Cl>>[CH3:9][CH:8]([CH3:10])[c:7]1[cH:11][cH:12][c:13]([cH:14][cH:15]1)[C:5](=[O:6])[CH2:4][CH2:3][CH2:2][Cl:1]"

    ReactionClassifier(rxn_smiles_with_atom_mapping, keep_mapping=True)


def test_get_template_smiles():
    rxn_smiles_with_atom_mapping = "[CH3:9][CH:8]([CH3:10])[c:7]1[cH:11][cH:12][cH:13][cH:14][cH:15]1.[CH2:3]([CH2:4][C:5](=[O:6])Cl)[CH2:2][Cl:1].[Al+3].[Cl-].[Cl-].[Cl-].C(Cl)Cl>>[CH3:9][CH:8]([CH3:10])[c:7]1[cH:11][cH:12][c:13]([cH:14][cH:15]1)[C:5](=[O:6])[CH2:4][CH2:3][CH2:2][Cl:1]"

    rxn_classifier = ReactionClassifier(rxn_smiles_with_atom_mapping, keep_mapping=True)

    assert rxn_classifier.get_template_smiles() == "c1ccccc1.CC(=O)Cl>>CC(=O)c1ccccc1"
