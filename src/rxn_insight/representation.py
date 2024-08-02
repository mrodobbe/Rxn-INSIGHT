from typing import Any

import numpy as np
import numpy.typing as npt
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdchem import Mol


def get_morgan_fingerprint(mol: Mol) -> npt.NDArray[Any]:
    """Get the ECFP4 fingerprint of a molecule.
    :param mol: RDKit Mol object
    :return: NumPy array
    """
    fp = AllChem.GetMorganFingerprintAsBitVect(
        mol, useChirality=True, radius=2, nBits=1024
    )
    morgan_fp = np.array(fp)
    return morgan_fp


def morgan_reaction_fingerprint(rxn: str) -> npt.NDArray[Any]:
    """Obtain the Morgan-based fingerprint of a reaction à la Schneider: https://doi.org/10.1021/ci5006614
    :param rxn: Reaction SMILES
    :return: NumPy array
    """
    reactants, products = rxn.split(">>")
    reactants_list = reactants.split(".")
    products_list = products.split(".")
    reactant_molecules = [Chem.AddHs(Chem.MolFromSmiles(r)) for r in reactants_list]
    product_molecules = [Chem.AddHs(Chem.MolFromSmiles(p)) for p in products_list]

    reactant_fp = tuple([get_morgan_fingerprint(mol) for mol in reactant_molecules])
    product_fp = tuple([get_morgan_fingerprint(mol) for mol in product_molecules])

    r_fp = np.sum(reactant_fp, axis=0)
    p_fp = np.sum(product_fp, axis=0)

    fp = p_fp - r_fp  # Difference fingerprint

    return fp


if __name__ == "__main__":
    rxn_smiles = "c1ccccc1.CC(=O)Cl>>CC(=O)c1ccccc1"
    difference_fingerprint = morgan_reaction_fingerprint(rxn_smiles)
    print(difference_fingerprint)
