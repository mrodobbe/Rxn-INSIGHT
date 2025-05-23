import hashlib
from typing import Any, Optional

import numpy as np
import numpy.typing as npt
import pandas as pd
from rdchiral.template_extractor import (
    MAXIMUM_NUMBER_UNMAPPED_PRODUCT_ATOMS,
    USE_STEREOCHEMISTRY,
    VERBOSE,
    canonicalize_transform,
    expand_changed_atom_tags,
    get_changed_atoms,
    get_fragments_for_changed_atoms,
    mols_from_smiles_list,
    replace_deuterated,
)
from rdkit import Chem
from rdkit.Chem import AllChem, MACCSkeys, rdChemReactions
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.rdchem import Atom, BondType, Mol
from rdkit.Chem.Scaffolds.MurckoScaffold import GetScaffoldForMol
from rxnmapper import RXNMapper
from scipy.spatial.distance import braycurtis, canberra, chebyshev, \
    cityblock, correlation, cosine, euclidean, minkowski, \
    dice, hamming, jaccard, kulczynski1, rogerstanimoto, \
    russellrao, sokalmichener, sokalsneath, yule

pd.options.mode.chained_assignment = None


def remove_atom_mapping(rxn: str, smarts: bool = False) -> str:
    """This function removes the mapping from mapped Reaction SMILES.
    :param smarts: SMIRKS instead of Reaction SMILES
    :param rxn: Reaction SMILES with mapping
    :return: Reaction SMILES without mapping
    """
    rxn = rxn.split(" |f")[0]  # Avoids errors with new style Reaction SMILES in USPTO
    reactants = rxn.split(">>")[0].split(
        "."
    )  # Make lists with all individual reactants
    products = rxn.split(">>")[1].split(".")  # Make lists with all individual products

    reactants_wo = []
    products_wo = []

    for r in reactants:
        if smarts:
            mol = Chem.MolFromSmarts(r)
        else:
            mol = Chem.MolFromSmiles(r)
        try:
            [a.SetAtomMapNum(0) for a in mol.GetAtoms()]
        except AttributeError:
            print(f"WARNING! {r} could not be parsed!")
            return ""  # Some molecules might disappear
        if smarts:
            s = Chem.MolToSmarts(mol)
        else:
            s = Chem.MolToSmiles(mol)
        reactants_wo.append(s)

    for p in products:
        if smarts:
            mol = Chem.MolFromSmarts(p)
        else:
            mol = Chem.MolFromSmiles(p)
        try:
            [a.SetAtomMapNum(0) for a in mol.GetAtoms()]
        except AttributeError:
            print(f"WARNING! {p} could not be parsed!")
            return ""  # Some molecules might disappear
        if smarts:
            s = Chem.MolToSmarts(mol)
        else:
            s = Chem.MolToSmiles(mol)
        products_wo.append(s)

    unmapped_rxn = f"{'.'.join(reactants_wo)}>>{'.'.join(products_wo)}"  # Combining unmapped reactants and products

    return unmapped_rxn


def get_atom_mapping(rxn: str, rxn_mapper: Optional[RXNMapper] = None) -> str:
    """This function maps reactants and products using RXNMapper (https://doi.org/10.1126/sciadv.abe4166)
    :param rxn_mapper: RXNMapper object
    :param rxn: Reaction SMILES without atom mapping
    :return: Reaction SMILES with atom mapping
    """
    if rxn_mapper is None:
        rxn_mapper = RXNMapper()

    try:
        results = rxn_mapper.get_attention_guided_atom_maps([rxn])[0]["mapped_rxn"]
        mapped_rxn: str = results
    except Exception as e:
        print(e)
        print("WARNING! Reaction could not be mapped! Returning unmapped reaction.")
        mapped_rxn = rxn

    return mapped_rxn


def sanitize_mapped_reaction(rxn: str) -> tuple[str, str, list[str]]:
    """Remove reactants that are unmapped from the reactants.
    :param rxn: Reaction SMILES with atom mapping
    :return: Mapped and unmapped reaction SMILES without reagents.
    """
    reactants = rxn.split(">>")[0].split(".")
    products = rxn.split(">>")[1].split(".")
    extra_agents = []

    mapped_reactants = []
    mapped_products = []
    for reactant in reactants:
        if ":" in reactant:
            if len(reactant) < 25:
                atoms = reactant.split("]")
                in_product = True
                for atom in atoms:
                    atom_map = atom.split(":")[-1]
                    if f":{atom_map}" not in rxn.split(">>")[1]:
                        in_product = False
                if in_product:
                    mapped_reactants.append(reactant)
                else:
                    extra_agents.append(reactant)
            else:
                mapped_reactants.append(reactant)
        else:
            extra_agents.append(reactant)

    for product in products:
        if ":" in product:
            mapped_products.append(product)
        else:
            extra_agents.append(product)

    # mapped_reactants = [reactant for reactant in reactants if ":" in reactant]
    # mapped_products = [product for product in products if ":" in product]

    for reactant in mapped_reactants:
        if reactant in mapped_products:
            mapped_reactants.remove(reactant)
            mapped_products.remove(reactant)
            extra_agents.append(reactant)

    sanitized_mapped_reaction = (
        ".".join(mapped_reactants) + ">>" + ".".join(mapped_products)
    )
    sanitized_unmapped_reaction = remove_atom_mapping(sanitized_mapped_reaction)

    return sanitized_mapped_reaction, sanitized_unmapped_reaction, extra_agents


def extract_from_reaction(
    reaction: dict[str, str | int], radius_reactants: int = 2, radius_products: int = 1
) -> dict[str, str | int]:
    """Extract the reaction template from mapped reaction SMILES. Code adapted
    from https://doi.org/10.1021/acs.jcim.9b00286.
    :param reaction: Dictionary with keys 'reactants' and 'products'.
    :param radius_reactants: Radius of atoms around the reaction center in the reactants
    :param radius_products: Radius of atoms around the reaction center in the products
    :return: dictionary with template information
    """
    reactants = mols_from_smiles_list(
        replace_deuterated(reaction["reactants"]).split(".")
    )
    products = mols_from_smiles_list(
        replace_deuterated(reaction["products"]).split(".")
    )

    # if rdkit cant understand molecule, return
    if None in reactants:
        return {"reaction_id": reaction["_id"]}
    if None in products:
        return {"reaction_id": reaction["_id"]}

    # try to sanitize molecules
    try:
        for i in range(len(reactants)):
            reactants[i] = AllChem.RemoveHs(reactants[i])  # *might* not be safe
        for i in range(len(products)):
            products[i] = AllChem.RemoveHs(products[i])  # *might* not be safe
        [Chem.SanitizeMol(mol) for mol in reactants + products]  # redundant w/ RemoveHs
        [mol.UpdatePropertyCache() for mol in reactants + products]
    except Exception as e:
        # can't sanitize -> skip
        print(e)
        print("Could not load SMILES or sanitize")
        print("ID: {}".format(reaction["_id"]))
        return {"reaction_id": reaction["_id"]}

    are_unmapped_product_atoms = False
    extra_reactant_fragment = ""
    for product in products:
        prod_atoms = product.GetAtoms()
        if sum([a.HasProp("molAtomMapNumber") for a in prod_atoms]) < len(prod_atoms):
            if VERBOSE:
                print("Not all product atoms have atom mapping")
            if VERBOSE:
                print("ID: {}".format(reaction["_id"]))
            are_unmapped_product_atoms = True

    if are_unmapped_product_atoms:  # add fragment to template
        for product in products:
            prod_atoms = product.GetAtoms()
            # Get unmapped atoms
            unmapped_ids = [
                a.GetIdx() for a in prod_atoms if not a.HasProp("molAtomMapNumber")
            ]
            if len(unmapped_ids) > MAXIMUM_NUMBER_UNMAPPED_PRODUCT_ATOMS:
                # Skip this example - too many unmapped product atoms!
                return {"reaction_id": reaction["_id"]}
            # Define new atom symbols for fragment with atom maps, generalizing fully
            atom_symbols = [f"[{a.GetSymbol()}]" for a in prod_atoms]
            # And bond symbols...
            bond_symbols = ["~" for _ in product.GetBonds()]
            if unmapped_ids:
                extra_reactant_fragment += (
                    AllChem.MolFragmentToSmiles(
                        product,
                        unmapped_ids,
                        allHsExplicit=False,
                        isomericSmiles=USE_STEREOCHEMISTRY,
                        atomSymbols=atom_symbols,
                        bondSymbols=bond_symbols,
                    )
                    + "."
                )
        if extra_reactant_fragment:
            extra_reactant_fragment = extra_reactant_fragment[:-1]
            if VERBOSE:
                print(f"    extra reactant fragment: {extra_reactant_fragment}")

        # Consolidate repeated fragments (stoichiometry)
        extra_reactant_fragment = ".".join(
            sorted(list(set(extra_reactant_fragment.split("."))))
        )

    if None in reactants + products:
        print("Could not parse all molecules in reaction, skipping")
        print("ID: {}".format(reaction["_id"]))
        return {"reaction_id": reaction["_id"]}

    # Calculate changed atoms
    changed_atoms, changed_atom_tags, err = get_changed_atoms(reactants, products)
    if err:
        if VERBOSE:
            print("Could not get changed atoms")
            print("ID: {}".format(reaction["_id"]))
        return {"reaction_id": reaction["_id"]}
    if not changed_atom_tags:
        if VERBOSE:
            print("No atoms changed?")
            print("ID: {}".format(reaction["_id"]))
        # print('Reaction SMILES: {}'.format(example_doc['RXN_SMILES']))
        return {"reaction_id": reaction["_id"]}

    try:
        # Get fragments for reactants
        reactant_fragments, intra_only, dimer_only = get_fragments_for_changed_atoms(
            reactants,
            changed_atom_tags,
            radius=radius_reactants,
            expansion=[],
            category="reactants",
        )
        # Get fragments for products
        # (WITHOUT matching groups but WITH the addition of reactant fragments)
        product_fragments, _, _ = get_fragments_for_changed_atoms(
            products,
            changed_atom_tags,
            radius=radius_products,
            expansion=expand_changed_atom_tags(changed_atom_tags, reactant_fragments),
            category="products",
        )
    except ValueError as e:
        if VERBOSE:
            print(e)
            print(reaction["_id"])
        return {"reaction_id": reaction["_id"]}

    # Put together and canonicalize (as good as possible)
    rxn_string = f"{reactant_fragments}>>{product_fragments}"
    rxn_canonical = canonicalize_transform(rxn_string)
    # Change from inter-molecular to intra-molecular
    rxn_canonical_split = rxn_canonical.split(">>")
    rxn_canonical = (
        rxn_canonical_split[0][1:-1].replace(").(", ".")
        + ">>"
        + rxn_canonical_split[1][1:-1].replace(").(", ".")
    )

    reactants_string = rxn_canonical.split(">>")[0]
    products_string = rxn_canonical.split(">>")[1]

    retro_canonical = products_string + ">>" + reactants_string

    # Load into RDKit
    rxn = AllChem.ReactionFromSmarts(retro_canonical)
    if rxn.Validate()[1] != 0:
        print("Could not validate reaction successfully")
        print("ID: {}".format(reaction["_id"]))
        print(f"retro_canonical: {retro_canonical}")
        # if VERBOSE:
        #     raw_input('Pausing...')
        return {"reaction_id": reaction["_id"]}

    template = {
        "products": products_string,
        "reactants": reactants_string,
        "reaction_smarts": retro_canonical,
        "intra_only": intra_only,
        "dimer_only": dimer_only,
        "reaction_id": reaction["_id"],
        "necessary_reagent": extra_reactant_fragment,
    }

    return template


def get_reaction_template(
    reaction: str, radius_reactants: int = 2, radius_products: int = 2
) -> str | None:
    """Get the reaction template from a mapped reaction.
    :param reaction: Mapped Reaction SMILES
    :param radius_reactants: Radius of atoms around the reaction center in the reactants
    :param radius_products: Radius of atoms around the reaction center in the products
    :return: Reaction template in forward direction
    """
    mapped_rxn: dict[str, str | int] = {
        "reactants": reaction.split(">>")[0],
        "products": reaction.split(">>")[1],
        "_id": 0,
    }

    result = extract_from_reaction(
        mapped_rxn, radius_reactants=radius_reactants, radius_products=radius_products
    )
    if "reactants" in result and "products" in result:
        template = f'{result["reactants"]}>>{result["products"]}'
        return template
    else:
        return None


def check_rings(atom: Atom, mol: Mol, match: list[int]) -> tuple[bool, list[int]]:
    if not atom.IsInRing():
        if len(list(atom.GetNeighbors())) + atom.GetNumExplicitHs() == 1:
            try:
                nb = atom.GetNeighbors()[0]
            except IndexError:
                return False, match
            if atom.GetAtomicNum() in [9, 17, 35, 53, 3, 11, 12, 29, 30, 34, 47, 50]:
                return False, match
            elif nb.IsInRing() and nb.GetIdx() in match:
                return True, match
            else:
                return False, match
        elif atom.GetAtomicNum() == 6 and (
            len(list(atom.GetNeighbors())) + atom.GetNumExplicitHs() != 1
        ):
            for nb in atom.GetNeighbors():
                if nb.IsInRing() and nb.GetAtomicNum() == 6 and nb.GetIdx() in match:
                    b = mol.GetBondBetweenAtoms(atom.GetIdx(), nb.GetIdx())
                    bo = b.GetBondTypeAsDouble()
                    if bo == 2:
                        return True, match
                    else:
                        return False, match
                else:
                    continue
            return False, match
        else:
            return False, match
    else:
        checker = False
        for nb in atom.GetNeighbors():
            if nb.GetIdx() in match:
                for ring in mol.GetRingInfo().AtomRings():
                    if nb.GetIdx() in ring and atom.GetIdx() in ring:
                        checker = True
                        match.append(atom.GetIdx())
                        break
    return checker, match


def atom_remover(mol: Mol, matches: list[list[int]]) -> Mol:
    if not matches:
        return Chem.Mol(mol)

    # We'll only process the first match (which is what the original code did)
    match = list(matches[0])

    # Create a new molecule
    res = Chem.RWMol(mol)
    res.BeginBatchEdit()

    # Collect ring information
    ring_info = mol.GetRingInfo()
    ring_atoms = set()
    for ring in ring_info.AtomRings():
        ring_atoms.update(ring)

    # Check if this is a ring-related operation
    is_ring_operation = all(idx in ring_atoms for idx in match)

    # If this is a ring operation, preserve ONLY the exact match atoms
    if is_ring_operation and len(match) >= 5:  # Assuming rings of size 5+
        # Remove atoms that are not in the match (strictly enforce ring atoms only)
        for aid in reversed(range(res.GetNumAtoms())):
            if aid not in match:
                res.RemoveAtom(aid)
    else:
        # For reaction templates, we need to preserve complete ring systems
        atoms_to_keep = set(match)

        # Add atoms that are part of the same ring as any atom in the match
        for ring in ring_info.AtomRings():
            if any(atom_idx in match for atom_idx in ring):
                # If yes, add all atoms in that ring to atoms_to_keep
                atoms_to_keep.update(ring)

        # Now remove atoms not in our keep list
        for aid in reversed(range(res.GetNumAtoms())):
            if aid not in atoms_to_keep:
                res.RemoveAtom(aid)

    res.CommitBatchEdit()

    try:
        Chem.SanitizeMol(res)
    except Exception:
        pass

    # For reaction templates, try to ensure proper ring representation
    if not is_ring_operation:
        try:
            smiles = Chem.MolToSmiles(res)
            new_mol = Chem.MolFromSmiles(smiles)
            if new_mol is not None:
                return new_mol
        except:
            pass

    return res


def draw_chemical_reaction(
    smiles: str, highlight_by_reactant: bool = False, font_scale: float = 1.5
) -> str:
    rxn = rdChemReactions.ReactionFromSmarts(smiles, useSmiles=True)
    trxn = rdChemReactions.ChemicalReaction(rxn)
    # move atom maps to be annotations:
    #     for m in trxn.GetReactants():
    #         moveAtomMapsToNotes(m)
    #     for m in trxn.GetProducts():
    #         moveAtomMapsToNotes(m)
    d2d = rdMolDraw2D.MolDraw2DSVG(800, 300)
    d2d.drawOptions().annotationFontScale = font_scale
    d2d.DrawReaction(trxn, highlightByReactant=highlight_by_reactant)

    d2d.FinishDrawing()
    drawing_text: str = d2d.GetDrawingText()
    return drawing_text


def move_atom_maps_to_notes(m: Mol) -> None:
    for at in m.GetAtoms():
        if at.GetAtomMapNum():
            at.SetProp("atomNote", str(at.GetAtomMapNum()))


def curate_smirks(df: pd.DataFrame) -> pd.DataFrame:
    """Make the SMIRKS database fit to the required format.
    :param df: Pandas DataFrame
    :return: Curated SMIRKS database
    """
    df["nreact"] = 0
    df["nproduct"] = 0
    for i in df.index:
        reactants = df["smirks"][i].split(">>")[0]
        df.loc[i, "nreact"] = len(reactants.split("."))
        products = df["smirks"][i].split(">>")[1]
        df.loc[i, "nproduct"] = len(products.split("."))
    return df


def get_map_index(mol: Mol) -> dict[int, int]:
    map_dict = dict()

    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        mapping = atom.GetAtomMapNum()
        map_dict[mapping] = idx

    return map_dict


def get_ring_systems(mol: Mol, include_spiro: bool = False) -> list[list[int]]:
    """Code taken from https://gist.github.com/greglandrum/de1751a42b3cae54011041dd67ae7415
    :param mol: RDKit Mol object
    :param include_spiro:
    :return: List with atoms that make up the ring systems
    """
    ri = mol.GetRingInfo()
    systems: list[list[int]] = []
    for ring in ri.AtomRings():
        ring_ats = set(ring)
        n_systems = []
        for system in systems:
            n_in_common = len(ring_ats.intersection(system))
            if n_in_common and (include_spiro or n_in_common > 1):
                ring_ats = ring_ats.union(system)
            else:
                n_systems.append(system)
        n_systems.append(list(ring_ats))
        systems = n_systems
    return systems


def remove_molecule_mapping(ring: Mol) -> str:
    ring = Chem.AddHs(ring)
    num_radicals = 0
    for atom in ring.GetAtoms():
        atom.SetAtomMapNum(0)
        num_radicals += atom.GetNumRadicalElectrons()
    try:
        ring = Chem.RemoveHs(ring)
    except:
        pass
    smiles: str = Chem.MolToSmiles(ring)
    if num_radicals > 0:
        smiles = smiles.replace("[CH]", "C")
        smiles = smiles.replace("[cH]", "c")
        smiles = smiles.replace("[NH]", "N")
        smiles = smiles.replace("[nH]", "n")
        smiles = smiles.replace("[OH]", "O")
        smiles = smiles.replace("[oH]", "o")
    return smiles


def sanitize_ring(mol: Mol) -> str:
    added_h = False
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
        nbs = atom.GetNeighbors()
        bos = []
        other_nbs = False
        num_nbs = 0
        ar_n = False
        for nb in nbs:
            num_nbs += 1
            b = mol.GetBondBetweenAtoms(atom.GetIdx(), nb.GetIdx())
            bos.append(b.GetBondTypeAsDouble())
            if b.GetBondTypeAsDouble() == 1.5 and atom.GetAtomicNum() == 7:
                ar_n = True
            if nb.GetAtomicNum() != 7:
                other_nbs = True
            if atom.GetIsAromatic() and nb.GetIsAromatic():
                b.SetBondType(BondType.AROMATIC)
        if atom.GetAtomicNum() == 7 and atom.GetIsAromatic():
            atom.SetNumExplicitHs(0)
        if (
            atom.GetAtomicNum() == 7
            and (sum(bos) != 3 or ar_n)
            and other_nbs
            and not added_h
            and num_nbs != 3
            and (not atom.GetIsAromatic() or not atom.IsInRingSize(6))
        ):
            atom.SetNumExplicitHs(1)
            added_h = True
        elif atom.GetAtomicNum() == 6 and atom.GetNumRadicalElectrons() != 0:
            explicit_hs = atom.GetNumExplicitHs()
            atom.SetNumExplicitHs(explicit_hs + atom.GetNumRadicalElectrons())
            atom.SetNumRadicalElectrons(0)
        if atom.GetAtomicNum() == 7 and atom.GetNumRadicalElectrons() != 0:
            atom.SetNumExplicitHs(1)

    # try:
    #     Chem.SanitizeMol(mol)
    # except Exception as e:
    #     pass

    new_mol = Chem.MolFromSmiles(Chem.MolToSmiles(mol))
    if new_mol is not None:
        smiles: str = Chem.MolToSmiles(mol)
        return smiles
    else:
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 7 and atom.GetNumExplicitHs() == 1:
                atom.SetNumExplicitHs(0)
                new_mol = Chem.MolFromSmiles(Chem.MolToSmiles(mol))
                if new_mol is not None:
                    smiles_explicit_H_fix: str = Chem.MolToSmiles(mol)
                    return smiles_explicit_H_fix
                else:
                    continue
            else:
                continue
        old_smiles: str = Chem.MolToSmiles(mol)
        return old_smiles


def get_scaffold(mol: Mol) -> str | None:
    """Get the Murcko scaffold of a molecule
    :param mol: RDKit Mol object
    :return: SMILES string
    """
    [a.SetAtomMapNum(0) for a in mol.GetAtoms()]
    scaffold = GetScaffoldForMol(mol)
    smi: Optional[str] = Chem.MolToSmiles(scaffold)

    return smi


def tag_reaction(rxn_info: dict[str, list[str] | str]) -> str:
    tag = f"{rxn_info['CLASS']} "
    fg_r = sorted(list(rxn_info["FG_REACTANTS"]))
    fg_p = sorted(list(rxn_info["FG_PRODUCTS"]))
    rings_r = sorted(list(rxn_info["PARTICIPATING_RINGS_REACTANTS"]))
    rings_p = sorted(list(rxn_info["PARTICIPATING_RINGS_PRODUCTS"]))
    tag += " ".join(fg_r) + " "
    tag += " ".join(fg_p) + " "
    tag += " ".join(rings_r) + " "
    tag += " ".join(rings_p)
    tag_bytes: bytes = tag.encode("UTF-8")
    hashtag = hashlib.sha256(tag_bytes).hexdigest()

    return str(hashtag)


def morgan_fp(mol: Mol) -> npt.NDArray[Any]:
    fp1 = AllChem.GetMorganFingerprintAsBitVect(
        mol, useChirality=True, radius=2, nBits=1024
    )
    vec1 = np.array(fp1)
    return vec1


def maccs_fp(mol: Chem.rdchem.Mol) -> npt.NDArray[Any]:
    return np.array(MACCSkeys.GenMACCSKeys(mol))


def get_fp(rxn: str, fp: str = "MACCS", concatenate: bool = True) -> npt.NDArray[Any]:
    reactant_str, product_str = rxn.split(">>")
    reactants = reactant_str.split(".")
    products = product_str.split(".")
    reactant_mols = [Chem.MolFromSmiles(reactant) for reactant in reactants]
    product_mols = [Chem.MolFromSmiles(product) for product in products]

    if fp.lower() == "maccs":
        reactant_fp = np.sum(np.array([maccs_fp(mol) for mol in reactant_mols]), axis=0)
        product_fp = np.sum(np.array([maccs_fp(mol) for mol in product_mols]), axis=0)
    elif fp.lower() == "morgan":
        reactant_fp = np.sum(
            np.array([morgan_fp(mol) for mol in reactant_mols]), axis=0
        )
        product_fp = np.sum(np.array([morgan_fp(mol) for mol in product_mols]), axis=0)
    else:
        raise KeyError(
            f"Fingerprint {fp} is not yet supported. Choose between MACCS and Morgan"
        )

    if concatenate:
        rxn_fp = np.concatenate((reactant_fp, product_fp))
    else:
        rxn_fp = np.sum((reactant_fp, product_fp), axis=0)

    return rxn_fp


def make_rdkit_fp(rxn: str, fp: str = "MACCS", concatenate: bool = True) -> str:
    nfp = get_fp(rxn, fp, concatenate)
    if len(np.where(nfp > 9)[0]) > 0:
        raise ValueError
    bfp = "".join(nfp.astype(str))
    return bfp


def get_similarity(v1: npt.NDArray[Any], v2: npt.NDArray[Any], metric: str = "jaccard") -> float:

    """
    Calculate the similarity between two fingerprints using a specified metric.

    Supported metrics include:
    - Binary metrics: `jaccard`, `dice`, `kulczynski1`, `rogerstanimoto`, `russellrao`, `sokalmichener`, `sokalsneath`, `yule`.
    - Distance-based metrics: `braycurtis`, `canberra`, `chebyshev`, `manhattan`, `correlation`, `cosine`, `euclidean`, `minkowski`.

    Args:
        v1 (npt.NDArray[Any]): Reference fingerprint as a NumPy array.
        v2 (npt.NDArray[Any]): Fingerprint to compare as a NumPy array.
        metric (str): Metric to calculate the similarity. Default is `"jaccard"`.

    Returns:
        float: Calculated similarity score.

    Raises:
        ValueError: If an unsupported metric is provided.

    Example:
        >>> import numpy as np
        >>> v1 = np.array([1, 0, 1, 1])
        >>> v2 = np.array([1, 1, 1, 0])
        >>> get_similarity(v1, v2, metric="jaccard")
        0.6666666666666667
    """

    metric = metric.lower()

    if metric == "jaccard":
        similarity = calculate_jaccard_similarity(v1, v2)
    elif metric == "dice":
        similarity = calculate_dice_similarity(v1, v2)
    elif metric == "kulczynski1":
        similarity = calculate_kulczynksi1_similarity(v1, v2)
    elif metric == "rogerstanimoto":
        similarity = calculate_rogerstanimoto_similarity(v1, v2)
    elif metric == "russellrao":
        similarity = calculate_russellrao_similarity(v1, v2)
    elif metric == "sokalmichener":
        similarity = calculate_sokalmichener_similarity(v1, v2)
    elif metric == "sokalsneath":
        similarity = calculate_sokalsneath_similarity(v1, v2)
    elif metric == "yule":
        similarity = calculate_yule_similarity(v1, v2)
    elif metric == "braycurtis":
        similarity = calculate_braycurtis_similarity(v1, v2)
    elif metric == "canberra":
        similarity = calculate_canberra_similarity(v1, v2)
    elif metric == "chebyshev":
        similarity = calculate_chebyshev_similarity(v1, v2)
    elif metric == "manhattan":
        similarity = calculate_manhattan_similarity(v1, v2)
    elif metric == "correlation":
        similarity = calculate_correlation_similarity(v1, v2)
    elif metric == "cosine":
        similarity = calculate_cosine_similarity(v1, v2)
    elif metric == "euclidean":
        similarity = calculate_euclidean_similarity(v1, v2)
    elif metric == "minkowski":
        similarity = calculate_minkowski_similarity(v1, v2)
    else:
        raise ValueError(f"Unknown metric '{metric}'. Please choose a valid metric.")

    return similarity


def calculate_jaccard_similarity(v1: npt.NDArray[Any], v2: npt.NDArray[Any]) -> float:
    """
    Calculate the Jaccard similarity between two vectors.

    The Jaccard similarity measures the proportion of shared elements in two sets.
    It is commonly used for binary vectors, indicating how similar the sets are.

    Args:
        v1 (npt.NDArray[Any]): First vector.
        v2 (npt.NDArray[Any]): Second vector.

    Returns:
        float: Jaccard similarity (1 - Jaccard distance).
    """
    return 1 - jaccard(v1, v2)


def calculate_dice_similarity(v1: npt.NDArray[Any], v2: npt.NDArray[Any]) -> float:
    """
    Calculate the Dice similarity between two vectors.

    Dice similarity is used for binary vectors and measures overlap between sets.

    Args:
        v1 (npt.NDArray[Any]): First vector.
        v2 (npt.NDArray[Any]): Second vector.

    Returns:
        float: Dice similarity (1 - Dice distance).
    """
    return 1 - dice(v1, v2)


def calculate_kulczynksi1_similarity(v1: npt.NDArray[Any], v2: npt.NDArray[Any]) -> float:
    """
    Calculate the Kulczynski 1 similarity between two vectors.

    Kulczynski similarity measures the average overlap between two sets.

    Args:
        v1 (npt.NDArray[Any]): First vector.
        v2 (npt.NDArray[Any]): Second vector.

    Returns:
        float: Kulczynski 1 similarity.
    """
    return 1 - kulczynski1(v1, v2)


def calculate_rogerstanimoto_similarity(v1: npt.NDArray[Any], v2: npt.NDArray[Any]) -> float:
    """
    Calculate the Rogers-Tanimoto similarity between two vectors.

    This metric evaluates similarity for binary vectors based on mismatched proportions.

    Args:
        v1 (npt.NDArray[Any]): First vector.
        v2 (npt.NDArray[Any]): Second vector.

    Returns:
        float: Rogers-Tanimoto similarity.
    """
    return 1 - rogerstanimoto(v1, v2)


def calculate_russellrao_similarity(v1: npt.NDArray[Any], v2: npt.NDArray[Any]) -> float:
    """
    Calculate the Russell-Rao similarity between two vectors.

    Russell-Rao similarity evaluates the proportion of matching 1's in binary vectors.

    Args:
        v1 (npt.NDArray[Any]): First vector.
        v2 (npt.NDArray[Any]): Second vector.

    Returns:
        float: Russell-Rao similarity.
    """
    return 1 - russellrao(v1, v2)


def calculate_sokalmichener_similarity(v1: npt.NDArray[Any], v2: npt.NDArray[Any]) -> float:
    """
    Calculate the Sokal-Michener similarity between two vectors.

    This metric is used for binary data, emphasizing matching 1's and 0's equally.

    Args:
        v1 (npt.NDArray[Any]): First vector.
        v2 (npt.NDArray[Any]): Second vector.

    Returns:
        float: Sokal-Michener similarity.
    """
    return 1 - sokalmichener(v1, v2)


def calculate_sokalsneath_similarity(v1: npt.NDArray[Any], v2: npt.NDArray[Any]) -> float:
    """
    Calculate the Sokal-Sneath similarity between two vectors.

    Sokal-Sneath similarity emphasizes 1's in binary data.

    Args:
        v1 (npt.NDArray[Any]): First vector.
        v2 (npt.NDArray[Any]): Second vector.

    Returns:
        float: Sokal-Sneath similarity.
    """
    return 1 - sokalsneath(v1, v2)


def calculate_yule_similarity(v1: npt.NDArray[Any], v2: npt.NDArray[Any]) -> float:
    """
    Calculate the Yule similarity between two vectors.

    Args:
        v1 (npt.NDArray[Any]): First vector.
        v2 (npt.NDArray[Any]): Second vector.

    Returns:
        float: Yule similarity.
    """
    return 1 - kulczynski1(v1, v2)


def calculate_braycurtis_similarity(v1: npt.NDArray[Any], v2: npt.NDArray[Any]) -> float:
    """
    Calculate the Bray-Curtis similarity between two vectors.

    This metric measures dissimilarity as the proportion of the absolute differences.

    Args:
        v1 (npt.NDArray[Any]): First vector.
        v2 (npt.NDArray[Any]): Second vector.

    Returns:
        float: Bray-Curtis similarity.
    """
    return 1 - braycurtis(v1, v2)


def calculate_canberra_similarity(v1: npt.NDArray[Any], v2: npt.NDArray[Any]) -> float:
    """
    Calculate the Canberra similarity between two vectors.

    The Canberra similarity is based on the sum of absolute differences scaled by the sum of vector elements.
    It emphasizes differences in smaller values.

    Args:
        v1 (npt.NDArray[Any]): First vector.
        v2 (npt.NDArray[Any]): Second vector.

    Returns:
        float: Canberra similarity (1 - Canberra distance).
    """
    return 1 - canberra(v1, v2)


def calculate_chebyshev_similarity(v1: npt.NDArray[Any], v2: npt.NDArray[Any]) -> float:
    """
    Calculate the Chebyshev similarity between two vectors.

    The Chebyshev similarity is derived from the maximum absolute difference between elements in the two vectors.

    Args:
        v1 (npt.NDArray[Any]): First vector.
        v2 (npt.NDArray[Any]): Second vector.

    Returns:
        float: Chebyshev similarity (1 - Chebyshev distance).
    """
    return 1 - chebyshev(v1, v2)


def calculate_manhattan_similarity(v1: npt.NDArray[Any], v2: npt.NDArray[Any]) -> float:
    """
    Calculate the Manhattan similarity between two vectors.

    The Manhattan similarity (or city block similarity) is the sum of the absolute differences between elements in the vectors.

    Args:
        v1 (npt.NDArray[Any]): First vector.
        v2 (npt.NDArray[Any]): Second vector.

    Returns:
        float: Manhattan similarity (1 - Manhattan distance).
    """
    return 1 - cityblock(v1, v2)


def calculate_correlation_similarity(v1: npt.NDArray[Any], v2: npt.NDArray[Any]) -> float:
    """
    Calculate the Correlation similarity between two vectors.

    The Correlation similarity measures the linear relationship between vectors, normalized to account for scale differences.

    Args:
        v1 (npt.NDArray[Any]): First vector.
        v2 (npt.NDArray[Any]): Second vector.

    Returns:
        float: Correlation similarity (1 - Correlation distance).
    """
    return 1 - correlation(v1, v2)


def calculate_cosine_similarity(v1: npt.NDArray[Any], v2: npt.NDArray[Any]) -> float:
    """
    Calculate the Cosine similarity between two vectors.

    The Cosine similarity measures the cosine of the angle between two vectors, indicating their directional similarity.

    Args:
        v1 (npt.NDArray[Any]): First vector.
        v2 (npt.NDArray[Any]): Second vector.

    Returns:
        float: Cosine similarity (1 - Cosine distance).
    """
    return 1 - cosine(v1, v2)


def calculate_euclidean_similarity(v1: npt.NDArray[Any], v2: npt.NDArray[Any]) -> float:
    """
    Calculate the Euclidean similarity between two vectors.

    The Euclidean similarity is derived from the straight-line distance between two points in multidimensional space.

    Args:
        v1 (npt.NDArray[Any]): First vector.
        v2 (npt.NDArray[Any]): Second vector.

    Returns:
        float: Euclidean similarity (1 - Euclidean distance).
    """
    return 1 - euclidean(v1, v2)


def calculate_minkowski_similarity(v1: npt.NDArray[Any], v2: npt.NDArray[Any], p: float = 3.0) -> float:
    """
    Calculate the Minkowski similarity between two vectors.

    The Minkowski similarity generalizes the Manhattan and Euclidean similarities using a parameter `p` to determine the metric.

    Args:
        v1 (npt.NDArray[Any]): First vector.
        v2 (npt.NDArray[Any]): Second vector.
        p (float): The order of the Minkowski metric (default: 3.0).

    Returns:
        float: Minkowski similarity (1 - Minkowski distance).
    """
    return 1 - minkowski(v1, v2, p=p)


def get_solvent_ranking(df: pd.DataFrame) -> pd.DataFrame:
    solvent_dict: dict[str, list[str]] = {"NAME": [], "COUNT": []}
    solvents = df["SOLVENT"].tolist()
    unique_solvents = list(set(solvents))
    for solvent in unique_solvents:
        solvent_dict["NAME"].append(solvent)
        solvent_dict["COUNT"].append(solvents.count(solvent))
    return pd.DataFrame(solvent_dict).sort_values(by="COUNT", ascending=False)


def get_catalyst_ranking(df: pd.DataFrame) -> pd.DataFrame:
    catalyst_dict: dict[str, list[str]] = {"NAME": [], "COUNT": []}
    catalysts = df["CATALYST"].tolist()
    unique_catalysts = list(set(catalysts))
    for catalyst in unique_catalysts:
        catalyst_dict["NAME"].append(catalyst)
        catalyst_dict["COUNT"].append(catalysts.count(catalyst))
    return pd.DataFrame(catalyst_dict).sort_values(by="COUNT", ascending=False)


def get_reagent_ranking(df: pd.DataFrame) -> pd.DataFrame:
    reagent_dict: dict[str, list[str]] = {"NAME": [], "COUNT": []}
    reagents = df["REAGENT"].tolist()
    unique_reagents = list(set(reagents))
    for reagent in unique_reagents:
        reagent_dict["NAME"].append(reagent)
        reagent_dict["COUNT"].append(reagents.count(reagent))
    return pd.DataFrame(reagent_dict).sort_values(by="COUNT", ascending=False)
