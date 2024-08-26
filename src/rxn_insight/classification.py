"""Reaction classification module"""

import itertools
from typing import Any, Optional, Union

import numpy as np
import numpy.typing as npt
import pandas as pd
from rdchiral.template_extractor import get_strict_smarts_for_atom
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdchem import Mol

try:
    from rxnmapper import RXNMapper
except ImportError:
    pass

from rxn_insight.utils import (
    atom_remover,
    get_atom_mapping,
    get_map_index,
    get_reaction_template,
    get_ring_systems,
    remove_atom_mapping,
    sanitize_mapped_reaction,
    sanitize_ring,
    tag_reaction,
)


class ReactionClassifier:
    """This class handles operations related to chemical reaction classification."""
    def __init__(
        self,
        reaction: str,
        rxn_mapper: Optional[RXNMapper] = None,
        keep_mapping: bool = False,
    ):
        """Initializes the ReactionClassifier with the specified reaction and options.

        Args:
            reaction (str): The reaction SMILES string with or without atom mapping.
            rxn_mapper (Optional[RXNMapper]): An instance of RXNMapper for generating atom mappings.
            keep_mapping (bool): If True, keeps existing atom mappings; otherwise, generates new mappings.
        """
        # Check reaction SMILES is valid
        try:
            AllChem.ReactionFromSmarts(reaction)
        except ValueError as e:
            raise ValueError(f"Invalid reaction SMILES string. Error msg: {e}")

        if keep_mapping:
            self.mapped_reaction = reaction
            self.reaction = remove_atom_mapping(
                reaction
            )  # Remove atom mapping for consistency
        else:
            self.reaction = reaction
            self.mapped_reaction = get_atom_mapping(
                self.reaction, rxn_mapper=rxn_mapper
            )
        self.sanitized_mapped_reaction, self.sanitized_reaction, self.extra_agents = (
            sanitize_mapped_reaction(self.mapped_reaction)
        )
        self.template = get_reaction_template(
            self.sanitized_mapped_reaction, radius_reactants=1, radius_products=0
        )
        self.template_smiles = self.get_template_smiles()
        self.reactants, self.products = self.sanitized_mapped_reaction.split(">>")
        self.reactant_mols = tuple(
            [Chem.MolFromSmiles(mol) for mol in self.reactants.split(".")]
        )
        self.product_mols = tuple(
            [Chem.MolFromSmiles(mol) for mol in self.products.split(".")]
        )
        self.mol_reactant = Chem.MolFromSmiles(self.reactants)
        self.mol_product = Chem.MolFromSmiles(self.products)
        self.n_atoms_reactants = self.mol_reactant.GetNumAtoms()
        self.n_atoms_products = self.mol_product.GetNumAtoms()
        self.num_reactants = len(self.sanitized_reaction.split(">>")[0].split("."))
        self.num_products = len(self.sanitized_reaction.split(">>")[1].split("."))
        try:
            (
                self.atom_mapping_index,
                self.atoms_diagonal,
                self.mapping_diagonal,
                self.matrix_size,
            ) = self.get_atom_mapping_indices()
        except ValueError:
            raise ValueError(
                "This reaction cannot be parsed, because no transformation is detected. "
                "Possibly, this is a resolution but RxnInsights cannot yet handle enantiomers."
            )
        self.be_matrix_reactants = self.get_be_matrix(self.mol_reactant)
        self.be_matrix_products = self.get_be_matrix(self.mol_product)
        self.r_matrix = self.be_matrix_products - self.be_matrix_reactants
        (
            self.transformation_matrix,
            self.reaction_center_atoms,
            self.transformation_mapping,
        ) = self.sanitize_r_matrix()
        self.removed_metals = False
        self.removed_halogens = False
        (
            self.sanitized_transformation_matrix,
            self.sanitized_reaction_center,
            self.sanitized_transformation_mapping,
        ) = self.remove_metals_and_halogens()
        self.nos_reaction_center = self.check_nos()
        self.ring_change = self.ring_changing()
        self.product_map_dict = get_map_index(self.mol_product)
        self.reactant_map_dict = get_map_index(self.mol_reactant)
        self.reaction_center_idx = [
            self.product_map_dict[i]
            for i in self.sanitized_transformation_mapping
            if i in self.product_map_dict
        ]
        self.transformation_determinant = np.linalg.det(
            self.sanitized_transformation_matrix
        )
        # self.reaction_class = self.classify_reaction()
        # self.functional_groups_reactants = self.get_functional_group_smarts(self.mol_reactant,
        #                                                                     self.be_matrix_reactants,
        #                                                                     self.reactant_map_dict)
        # self.functional_groups_products = self.get_functional_group_smarts(self.mol_product,
        #                                                                    self.be_matrix_products,
        #                                                                    self.product_map_dict)

    def get_template_smiles(self) -> str | None:
        """Generates a reaction SMILES from the reaction SMARTS template.

        Returns:
            str | None: The reaction SMILES of the reaction template, or None if no template is generated.
        """
        extended_template = self.template
        if extended_template is None:
            return None
        reaction = self.sanitized_mapped_reaction

        reactants = reaction.split(">>")[0].split(".")
        products = reaction.split(">>")[1].split(".")

        reactants_template, products_template = extended_template.split(">>")
        reactants_template_list = reactants_template.split(".")
        products_template_list = products_template.split(".")

        subreactants = []
        subproducts = []

        for j in range(len(reactants)):
            reactant_molecule = Chem.MolFromSmiles(reactants[j])
            atoms_reactant = []
            for g in range(len(reactants_template_list)):
                reactant_template = Chem.MolFromSmarts(reactants_template_list[g])
                if len(reactant_molecule.GetSubstructMatch(reactant_template)) > 0:
                    atoms_reactant.append(
                        reactant_molecule.GetSubstructMatch(reactant_template)
                    )
                    break
            m = atom_remover(reactant_molecule, atoms_reactant)
            for atom in m.GetAtoms():
                explicit_hs = atom.GetNumExplicitHs()
                atom.SetNumExplicitHs(explicit_hs + atom.GetNumRadicalElectrons())
                atom.SetNumRadicalElectrons(0)
                atom.SetAtomMapNum(0)
            subreactants.append(m)

        for k in range(len(products)):
            #         print(products[k])
            product_molecule = Chem.MolFromSmiles(products[k])
            atoms_products = []
            for g in range(len(products_template_list)):
                product_template = Chem.MolFromSmarts(products_template_list[g])
                if len(product_molecule.GetSubstructMatch(product_template)) > 0:
                    atoms_products.append(
                        product_molecule.GetSubstructMatch(product_template)
                    )
                    break
            m = atom_remover(product_molecule, atoms_products)
            for atom in m.GetAtoms():
                explicit_hs = atom.GetNumExplicitHs()
                atom.SetNumExplicitHs(explicit_hs + atom.GetNumRadicalElectrons())
                atom.SetNumRadicalElectrons(0)
                atom.SetAtomMapNum(0)
            subproducts.append(m)

        rxn: str = (
            ".".join([Chem.MolToSmiles(r) for r in subreactants])
            + ">>"
            + Chem.MolToSmiles(subproducts[0])
        )

        return rxn

    def get_functional_group_smarts(
        self, molecule: Mol, matrix: npt.NDArray[Any], map_dict: dict[int, int]
    ) -> tuple[str, ...]:
        """Identifies and returns SMARTS strings for functional groups in the molecule based on the specified matrix and mapping.

        Args:
            molecule (Mol): The RDKit molecule object.
            matrix (npt.NDArray[Any]): A matrix representing chemical properties or structure.
            map_dict (dict[int, int]): Mapping of atom indices to their corresponding mapping numbers in the reaction.

        Returns:
            tuple[str, ...]: A tuple containing SMARTS strings of the identified functional groups.
        """
        maps = self.transformation_mapping
        matrix_indices = [self.atom_mapping_index[atom_map] for atom_map in maps]
        functional_groups = []
        bond_order_dict = {1.0: "-", 1.5: ":", 2.0: "=", 3.0: "#"}
        idx_map_dict = {v: k for k, v in self.atom_mapping_index.items()}
        visited_atoms = np.array([])
        matrix_indices_np = np.copy(np.array(matrix_indices))
        del_idx = 0
        for i in range(len(matrix_indices)):
            idx = matrix_indices[i]
            if self.atoms_diagonal[idx] == 6 and not np.all(
                self.atoms_diagonal[matrix_indices[i:]] == 0
            ):
                matrix_indices_np = np.delete(matrix_indices_np, del_idx)
                matrix_indices_np = np.append(matrix_indices_np, idx)
            else:
                del_idx += 1

        matrix_indices = matrix_indices_np

        for i in range(len(matrix_indices)):
            idx = matrix_indices[i]
            r = matrix[:, idx]
            vals = np.array(list(set(list(np.where(r != 0.0)[0]) + [idx])))
            vals_new = []
            atom_maps = []
            atom_indices = []
            for val in vals:
                if val in idx_map_dict:
                    mapping = idx_map_dict[val]
                    if mapping in map_dict:
                        vals_new.append(val)
                        atom_maps.append(mapping)
                        atom_indices.append(map_dict[mapping])
            vals = np.array(vals_new)
            if len(visited_atoms[np.in1d(visited_atoms, vals)]) > 0:
                continue
            elif len(vals) == 0:
                continue

            fg_matrix = matrix[vals][:, vals]
            try:
                matrix_id = np.where(np.array(vals) == idx)[0][0]
            except IndexError:
                continue
            bond_orders = fg_matrix[:, matrix_id]
            main_atom_id = atom_indices[matrix_id]
            main_atom = molecule.GetAtomWithIdx(main_atom_id)
            smarts_mapped = get_strict_smarts_for_atom(main_atom)
            if ":" in smarts_mapped:
                smarts = smarts_mapped.split(":")[0] + "]"
            else:
                smarts = smarts_mapped
            neighbors_to_go = len(bond_orders) - 1

            sort_bonds = (-bond_orders).argsort()
            bond_orders = bond_orders[sort_bonds]
            vals = vals[sort_bonds].astype(np.int32)
            atom_indices = np.array(atom_indices)[sort_bonds]

            for j in range(len(bond_orders)):
                if vals[j] == idx:
                    continue
                elif bond_orders[j] == 0:
                    neighbors_to_go -= 1
                    continue
                else:
                    atom = molecule.GetAtomWithIdx(int(atom_indices[j]))
                    bond = bond_order_dict[bond_orders[j]]
                    smarts_mapped = get_strict_smarts_for_atom(atom)
                    if ":" in smarts_mapped:
                        smarts_unmapped = smarts_mapped.split(":")[0] + "]"
                    else:
                        smarts_unmapped = smarts_mapped
                    if neighbors_to_go == 1:
                        smarts += f"{bond}{smarts_unmapped}"
                    else:
                        smarts += f"({bond}{smarts_unmapped})"
                        neighbors_to_go -= 1
            functional_groups.append(smarts)
            visited_atoms = np.append(visited_atoms, vals)

        return tuple(functional_groups)

    def get_functional_groups(
        self, mol: Mol, map_dict: dict[int, int], df: pd.DataFrame
    ) -> list[str]:
        """Extracts functional groups from the molecule using the specified mapping and reference DataFrame.

        Args:
            mol (Mol): The molecule from which to extract functional groups.
            map_dict (dict[int, int]): A dictionary mapping atom indices to mapping numbers.
            df (pd.DataFrame): DataFrame containing functional group definitions.

        Returns:
            list[str]: A list of names of identified functional groups.
        """
        maps = self.transformation_mapping
        atom_indices = np.array(
            [map_dict[atom_map] for atom_map in maps if atom_map in map_dict]
        )
        fg = []
        visited_atoms: list[list[int]] = []
        for i in df.index:
            if len(np.in1d(visited_atoms, atom_indices)) != 0:
                if len(visited_atoms[np.in1d(visited_atoms, atom_indices)]) == len(
                    atom_indices
                ):
                    break
            sm = mol.GetSubstructMatches(Chem.MolFromSmarts(df["pattern"][i]))
            if len(sm) == 0:
                continue
            else:
                for m in sm:
                    matched_atoms = np.array(m)
                    if len(matched_atoms[np.in1d(matched_atoms, atom_indices)]) > 0:
                        if len(np.in1d(visited_atoms, matched_atoms)) == 0:
                            fg.append(df["name"][i])
                            visited_atoms = np.unique(
                                np.append(visited_atoms, matched_atoms)
                            )
                        elif len(
                            visited_atoms[np.in1d(visited_atoms, matched_atoms)]
                        ) != len(matched_atoms):
                            fg.append(df["name"][i])
                            visited_atoms = np.unique(
                                np.append(visited_atoms, matched_atoms)
                            )
                        else:
                            continue
                    else:
                        continue
        return fg

    def get_ring_type(
        self, mol: Mol, map_dict: Optional[dict[int, int]] = None
    ) -> list[str]:
        """Determines the types of ring structures present in the molecule.

        Args:
            mol (Mol): The molecule to analyze.
            map_dict (Optional[dict[int, int]]): Mapping of atom indices to their mapping numbers, if available.

        Returns:
            list[str]: A list of ring types identified in the molecule.
        """
        try:
            rs = get_ring_systems(mol, include_spiro=True)
        except:
            return []

        if map_dict is not None:
            if len(rs) == 0:
                return []
            else:
                involved_rings = []
                maps = self.transformation_mapping
                atom_indices = np.array(
                    [map_dict[atom_map] for atom_map in maps if atom_map in map_dict]
                )
                for r in rs:
                    r = np.array(r)
                    if np.in1d(r, atom_indices).sum() > 0:
                        involved_rings.append(r)

            if len(involved_rings) == 0:
                return []
            else:
                rs = involved_rings

        found_rings = []
        if len(rs) > 0:
            for k in range(len(rs)):
                found_rings.append(sanitize_ring(atom_remover(mol, [rs[k]])))
            return found_rings
        else:
            return []

    def balance_reaction(self, fgr: list[str], fgp: list[str]) -> list[str]:
        """Balances the reaction based on functional groups present in reactants and products.

        Args:
            fgr (list[str]): Functional groups in reactants.
            fgp (list[str]): Functional groups in products.

        Returns:
            list[str]: A list of potential by-products or missing elements in the balanced reaction.
        """
        d = self.transformation_matrix.diagonal()
        mr = self.be_matrix_reactants
        mp = self.be_matrix_products
        lost_heavy = self.mol_reactant.GetNumAtoms() - self.mol_product.GetNumAtoms()
        if lost_heavy == 0:
            return []
        negative_values = np.where(d < 0)[0]
        metals = np.array([3, 5, 11, 12, 29, 30, 34, 47, 50])
        metal_indices = np.where(np.in1d(self.reaction_center_atoms, metals))[0]
        negative_values = np.unique(
            np.array(list(negative_values) + list(metal_indices))
        )
        if len(negative_values) == 0:
            return ["Other"]
        atom_indices = self.transformation_mapping[negative_values]
        matrix_indices = np.array([self.atom_mapping_index[a] for a in atom_indices])
        if len(matrix_indices) == 0:
            return ["Other"]
        lost_atoms = []
        for idx in matrix_indices:
            in_reactants = len(np.where(mr[:, idx] != 0)[0]) > 0
            in_products = len(np.where(mp[:, idx] != 0)[0]) > 0
            if in_reactants and not in_products:
                symbol = int(self.atoms_diagonal[idx])
                lost_atoms.append(Chem.GetPeriodicTable().GetElementSymbol(symbol))
        try:
            lost_all = (
                Chem.AddHs(self.mol_reactant).GetNumAtoms()
                - Chem.AddHs(self.mol_product).GetNumAtoms()
            )
            n_lost_hs = lost_all - lost_heavy
            lost_hs = ["H" for _ in range(n_lost_hs)]
        except:
            print("WARNING! Could not calculate lost hydrogens.")
            lost_hs = []
        all_lost = lost_atoms + lost_hs
        unique_lost = list(set(all_lost))
        lost_dict = {"H": 0}
        for atom in unique_lost:
            if atom in lost_dict:
                lost_dict[atom] += 1
            else:
                lost_dict[atom] = 1

        small_molecules = []

        if lost_heavy > 7 and "Tosylate" in fgr and "Tosylate" not in fgp:
            small_molecules.append("TsOH")
            if "S" in lost_dict:
                lost_dict["S"] -= 1
            if "O" in lost_dict:
                lost_dict["O"] -= 3
        elif lost_heavy > 7 and "Triflate" in fgr and "Triflate" not in fgp:
            small_molecules.append("TfOH")
            if "S" in lost_dict:
                lost_dict["S"] -= 1
            if "O" in lost_dict:
                lost_dict["O"] -= 3
            if "F" in lost_dict:
                lost_dict["F"] -= 3
        elif lost_heavy > 4 and "Mesylate" in fgr and "Mesylate" not in fgp:
            small_molecules.append("MsOH")
            if "S" in lost_dict:
                lost_dict["S"] -= 1
            if "O" in lost_dict:
                lost_dict["O"] -= 3

        for x in ["F", "Cl", "Br", "I"]:
            if x not in lost_dict:
                continue
            elif lost_dict["H"] > 0 and lost_dict[x] > 0:
                small_molecules.append(f"H{x}")
                lost_dict["H"] -= 1
                lost_dict[x] -= 1
            elif x == "Br" or x == "I" and lost_dict[x] > 0:
                small_molecules.append(f"{x}-")
                lost_dict[x] -= 1
                continue
        for x in ["O", "S"]:
            if x not in lost_dict:
                continue
            elif lost_dict["H"] > 1:
                small_molecules.append(f"H2{x}")
                lost_dict["H"] -= 2
                lost_dict[x] -= 1
            elif x == "O" and lost_dict["H"] == 1 and "S" not in lost_dict:
                small_molecules.append("HO-")
                lost_dict["H"] -= 1
        if "N" in lost_dict and lost_dict["H"] > 2:
            small_molecules.append(f"H2{x}")
            lost_dict["H"] -= 3
            lost_dict["N"] -= 1
            small_molecules.append("NH3")
        for x in [
            Chem.GetPeriodicTable().GetElementSymbol(int(symbol)) for symbol in metals
        ]:
            if x not in lost_dict:
                continue
            else:
                small_molecules.append(x)

        if lost_heavy > 0 and len(small_molecules) == 0:
            small_molecules = ["Other"]

        return small_molecules

    def get_reaction_center_info(self, df: pd.DataFrame) -> dict[str, Union[list[str], str, int]]:
        """Compiles detailed information about the reaction center from the reaction.

        Args:
            df (pd.DataFrame): DataFrame containing additional data required for analysis.

        Returns:
            dict[str, Union[list[str], str, int]]:
            A dictionary containing detailed information about the reaction center.
        """
        reaction_center: dict[str, Union[list[str], str, int]] = dict()
        reaction_center["REACTION"] = self.sanitized_reaction
        reaction_center["MAPPED_REACTION"] = self.sanitized_mapped_reaction
        reaction_center["N_REACTANTS"] = self.num_reactants
        reaction_center["N_PRODUCTS"] = self.num_products
        fg_reactants = self.get_functional_groups(
            self.mol_reactant, self.reactant_map_dict, df
        )
        reaction_center["FG_REACTANTS"] = fg_reactants
        fg_products = self.get_functional_groups(
            self.mol_product, self.product_map_dict, df
        )
        reaction_center["FG_PRODUCTS"] = fg_products
        reaction_center["PARTICIPATING_RINGS_REACTANTS"] = self.get_ring_type(
            self.mol_reactant, self.reactant_map_dict
        )
        reaction_center["PARTICIPATING_RINGS_PRODUCTS"] = self.get_ring_type(
            self.mol_product, self.product_map_dict
        )
        reaction_center["ALL_RINGS_PRODUCTS"] = self.get_ring_type(self.mol_product)
        reaction_center["BY-PRODUCTS"] = self.balance_reaction(
            fg_reactants, fg_products
        )
        reaction_center["CLASS"] = self.classify_reaction()
        reaction_center["TAG"] = tag_reaction(reaction_center)
        return reaction_center

    def get_atom_mapping_indices(
        self,
    ) -> tuple[dict[int, int], npt.NDArray[Any], npt.NDArray[Any], int]:
        """Generates a mapping from atom indices to their positions in the transformation matrix.

        Returns:
            tuple: Contains a dictionary for atom mapping to indices, arrays for atom numbers and mapping, and the matrix size.
        """

        """Make a dictionary that gives a unique index to all atoms in reactants and products.
        Necessary since reactions are not balanced.
        :return: Dictionary that links atom map and index. Size of BE-matrix
        """
        map_idx_dict = dict()
        atom_number_dict = dict()
        unmapped_reactant_atoms = []
        i = 0
        for atom in self.mol_reactant.GetAtoms():
            atom_idx = atom.GetIdx()
            atom_map = atom.GetAtomMapNum()
            if atom_map == 0:
                unmapped_reactant_atoms.append(atom_idx)
                continue
            else:
                map_idx_dict[atom_map] = i
                atom_number_dict[i] = atom.GetAtomicNum()
                i += 1

        unmapped_product_atoms = []
        for atom in self.mol_product.GetAtoms():
            atom_idx = atom.GetIdx()
            atom_map = atom.GetAtomMapNum()
            if atom_map == 0:
                unmapped_product_atoms.append(atom_idx)
                continue
            elif atom_map in map_idx_dict:
                continue
            else:
                map_idx_dict[atom_map] = i
                atom_number_dict[i] = atom.GetAtomicNum()
                i += 1

        atom_mapping_values = list(map_idx_dict.keys())
        max_map = max(atom_mapping_values) + 1

        for unmapped_atom in unmapped_reactant_atoms:
            atom = self.mol_reactant.GetAtomWithIdx(unmapped_atom)
            atom.SetAtomMapNum(max_map)
            map_idx_dict[max_map] = i
            atom_number_dict[i] = atom.GetAtomicNum()
            max_map += 1
            i += 1

        for unmapped_atom in unmapped_product_atoms:
            atom = self.mol_product.GetAtomWithIdx(unmapped_atom)
            atom.SetAtomMapNum(max_map)
            map_idx_dict[max_map] = i
            atom_number_dict[i] = atom.GetAtomicNum()
            max_map += 1
            i += 1

        matrix_size = len(map_idx_dict.keys())
        atoms_diagonal = np.array([atom_number_dict[j] for j in range(matrix_size)])
        idx_map_dict = {v: k for k, v in map_idx_dict.items()}
        map_diagonal = np.array([idx_map_dict[j] for j in range(matrix_size)])

        return map_idx_dict, atoms_diagonal, map_diagonal, matrix_size

    def get_be_matrix(self, molecule: Mol) -> npt.NDArray[Any]:
        """Calculates the bond-electron matrix for the given molecule.

        Args:
            molecule (Mol): The molecule for which to calculate the bond-electron matrix.

        Returns:
            npt.NDArray[Any]: A matrix representing the bond-electron relationships in the molecule.
        """
        matrix = np.zeros((self.matrix_size, self.matrix_size))
        for atom in molecule.GetAtoms():
            atom_map = atom.GetAtomMapNum()
            atom_idx = atom.GetIdx()
            idx = self.atom_mapping_index[atom_map]
            matrix[idx][idx] = (
                Chem.GetPeriodicTable().GetNOuterElecs(atom.GetAtomicNum())
                - atom.GetTotalValence()
            )

            for nb in atom.GetNeighbors():
                nb_map = nb.GetAtomMapNum()
                nb_idx = self.atom_mapping_index[nb_map]
                nb_atom_idx = nb.GetIdx()
                bond = molecule.GetBondBetweenAtoms(atom_idx, nb_atom_idx)
                if bond is not None:
                    matrix[idx][nb_idx] = bond.GetBondTypeAsDouble()
                    matrix[nb_idx][idx] = matrix[idx][nb_idx]

        return matrix

    def sanitize_r_matrix(self) -> tuple[npt.NDArray[Any], ...]:
        """Sanitizes the R-matrix by removing all-zero rows and columns.

        Returns:
            tuple[npt.NDArray[Any], ...]: A tuple containing the cleaned R-matrix and arrays for atom numbers and mappings.
        """
        clean_r_matrix = np.copy(self.r_matrix)
        atoms_to_remove = ~np.all(clean_r_matrix == 0.0, axis=1)
        clean_r_matrix = clean_r_matrix[atoms_to_remove][
            :, ~np.all(clean_r_matrix[atoms_to_remove] == 0.0, axis=0)
        ]
        sanitized_atoms_diagonal = self.atoms_diagonal[atoms_to_remove]
        sanitized_mapping_diagonal = self.mapping_diagonal[atoms_to_remove]

        return clean_r_matrix, sanitized_atoms_diagonal, sanitized_mapping_diagonal

    def remove_metals_and_halogens(self) -> tuple[npt.NDArray[Any], ...]:
        """Removes metal and halogen atoms from the R-matrix, as they are generally not useful for reaction classification.

        Returns:
            tuple[npt.NDArray[Any], ...]: A tuple containing the sanitized transformation matrix, reaction center atoms, and mappings.
        """
        to_delete = []
        ignore_atoms = [9, 17, 35, 53]
        metals = [3, 5, 11, 12, 29, 30, 34, 47, 50]
        for i in range(len(self.reaction_center_atoms)):
            if (
                self.reaction_center_atoms[i] in ignore_atoms
                and self.transformation_matrix.diagonal()[i] < 0
            ):
                if i not in to_delete:
                    self.removed_halogens = True
                    to_delete.append(i)
            elif (
                self.reaction_center_atoms[i] in metals
                and len(np.where(np.array(self.transformation_matrix[:, i]) < 0)[0]) > 0
            ):
                if i not in to_delete:
                    self.removed_metals = True
                    to_delete.append(i)

        to_delete = sorted(to_delete)

        clean_transformation_matrix = np.copy(self.transformation_matrix)
        reaction_center = np.copy(self.reaction_center_atoms)
        transformation_mapping = np.copy(self.transformation_mapping)
        for j in reversed(range(len(to_delete))):
            clean_transformation_matrix = np.delete(
                clean_transformation_matrix, (to_delete[j]), axis=0
            )
            clean_transformation_matrix = np.delete(
                clean_transformation_matrix, (to_delete[j]), axis=1
            )
            reaction_center = np.delete(reaction_center, (to_delete[j]))
            transformation_mapping = np.delete(transformation_mapping, (to_delete[j]))

        return clean_transformation_matrix, reaction_center, transformation_mapping

    def check_nos(self) -> bool:
        """Checks if nitrogen, oxygen, or sulfur atoms are involved in the reaction center.

        Returns:
            bool: True if N, O, or S atoms are involved in the reaction center; otherwise, False.
        """
        nos = False
        nos_atoms = [7, 8, 16]
        for i in range(len(self.sanitized_reaction_center)):
            atom = self.sanitized_reaction_center[i]
            if (
                atom in nos_atoms
                and len(np.where(self.sanitized_transformation_matrix[:, i] > 0)[0]) > 0
            ):
                nos = True
                break

        return nos

    def ring_changing(self) -> int:
        """Calculates the net change in the number of ring structures between reactants and products.

        Returns:
            int: The net change in the number of rings; positive for ring formation, negative for ring breaking.
        """
        reactants = Chem.AddHs(self.mol_reactant)
        products = Chem.AddHs(self.mol_product)

        ri_r: int = reactants.GetRingInfo().NumRings()
        ri_p: int = products.GetRingInfo().NumRings()

        ri_change = ri_p - ri_r

        return ri_change

    def is_fgi(self) -> bool:
        """Determines if the reaction involves a functional group interconversion (FGI).

        Returns:
            bool: True if the reaction is classified as a functional group interconversion, otherwise False.
        """
        m = self.transformation_matrix
        if len(self.sanitized_transformation_matrix) == 1 and m[0][0] == 0:
            return True
        d = m.diagonal()
        a = np.array(
            [
                self.atom_mapping_index[atom_map]
                for atom_map in self.transformation_mapping
            ]
        )
        if len(a) == 0:
            return False
        involved_types = self.atoms_diagonal[a]
        indices = np.where(involved_types == 6)[0]
        s_indices = np.where(involved_types == 16)[0]
        si_indices = np.where(involved_types == 14)[0]
        o_indices = np.where(involved_types == 8)[0]
        x_indices = np.array(
            list(np.where(involved_types == 9)[0])
            + list(np.where(involved_types == 17)[0])
            + list(np.where(involved_types == 35)[0])
            + list(np.where(involved_types == 53)[0])
        )
        for idx in indices:
            if sum(m[:, idx]) != 0:
                if len(np.where(d > 0)[0]) == len(np.where(d < 0)[0]) and len(
                    np.where(m[:, idx] > 0)[0]
                ) == len(np.where(m[:, idx] < 0)[0]):
                    return True
                otf = "[O;D2;+0]-[S;D4;+0](=[O;H0;D1;+0])(=[O;H0;D1;+0])-[CX4;D4;+0](-F)(-F)-F"
                ots = (
                    "[C;H3;D1;+0]-[c;H0;D3;+0]1:[c;H1;D2;+0]:[c;H1;D2;+0]:[c;H0;D3;+0]"
                    "(-[S;H0;D4;+0](=[O;H0;D1;+0])(=[O;H0;D1;+0])):[c;H1;D2;+0]:[c;H1;D2;+0]:1"
                )
                oms = "[O;H0;D2;+0]-[S;H0;D4;+0](=[O;H0;D1;+0])(=[O;H0;D1;+0])-[C;H3;D1;+0]"
                thiocyanate = "[*]-[S;H0;D2;+0]-[C;H0;D2;+0]#[N;H0;D1;+0]"
                isothiocyanate = "[*]-[N;H0;D2;+0]=[C;H0;D2;+0]=[S;H0;D1;+0]"
                isocyanide = "[#6]-[N;H0;D2;+1]#[C;H0;D1;-1]"
                isocyanate = "[#6]-[N;H0;D2;+0]=[C;H0;D2;+0]=[O;H0;D1;+0]"

                leaving = [
                    Chem.MolFromSmarts(otf),
                    Chem.MolFromSmarts(ots),
                    Chem.MolFromSmarts(oms),
                    Chem.MolFromSmarts(thiocyanate),
                    Chem.MolFromSmarts(isothiocyanate),
                    Chem.MolFromSmarts(isocyanate),
                    Chem.MolFromSmarts(isocyanide),
                ]
                for group in leaving:
                    if (
                        len(self.mol_reactant.GetSubstructMatches(group)) > 0
                        and len(self.mol_product.GetSubstructMatches(group)) == 0
                    ):
                        if self.is_cc_coupling() or self.is_heteroatom_alkylation():
                            return False
                        else:
                            return True
                    elif (
                        len(self.mol_reactant.GetSubstructMatches(group)) > 0
                        and len(self.mol_product.GetSubstructMatches(group)) > 0
                    ):
                        if (
                            len(
                                a[
                                    np.in1d(
                                        a, self.mol_product.GetSubstructMatches(group)
                                    )
                                ]
                            )
                            > 0
                        ):
                            return True
                        else:
                            continue
                    else:
                        continue
                if (
                    len(np.where(d > 0)[0]) == 0
                    and len(np.where(d < 0)[0]) > 0
                    and len(np.where(m[:, idx] < 0)[0]) == 0
                    and self.num_reactants == 1
                ):
                    return True
                elif (
                    len(np.where(d < 0)[0]) == 0
                    and len(np.where(d > 0)[0]) > 0
                    and len(np.where(m[:, idx] > 0)[0]) == 0
                    and self.num_reactants == 1
                    and not self.is_oxidation()
                ):
                    return True
                elif (
                    len(np.where(d > 0)[0]) == 0
                    and len(np.where(d < 0)[0]) > 0
                    and self.num_reactants == 2
                    and len(x_indices) > 0
                ):
                    for x_idx in x_indices:
                        if len(np.where(m[:, x_idx] < 0)[0]) == len(
                            np.where(m[:, x_idx] > 0)[0]
                        ):
                            return True
                else:
                    boronic = Chem.MolFromSmarts("[#6]-[BX3]")
                    nitro = Chem.MolFromSmarts(
                        "[$([NX3](=O)=O),$([NX3+](=O)[O-])][!#8]"
                    )
                    if (
                        len(self.mol_product.GetSubstructMatches(boronic)) > 0
                        and len(self.mol_reactant.GetSubstructMatches(boronic)) == 0
                    ):
                        return True
                    elif (
                        len(self.mol_product.GetSubstructMatches(nitro)) == 0
                        and 0
                        < len(self.mol_reactant.GetSubstructMatches(nitro))
                        < len(np.where(m > 0)[0])
                        and not self.is_reduction()
                    ):
                        return True
                    else:
                        return False
        if d.any():
            nitro = Chem.MolFromSmarts("[$([NX3](=O)=O),$([NX3+](=O)[O-])][!#8]")
            if len(self.mol_product.GetSubstructMatches(nitro)) == 0 and 0 < len(
                self.mol_reactant.GetSubstructMatches(nitro)
            ):
                if not self.is_reduction():
                    return True
                else:
                    return False
            if len(np.where(d > 0)[0]) == len(np.where(d < 0)[0]):
                if len(s_indices) > 0 and len(o_indices) > 0 and self.is_oxidation():
                    return False
                else:
                    return True
            elif len(np.where(d > 0)[0]) > 0 and len(np.where(d < 0)[0]) > 0:
                if len(s_indices) > 0 and len(o_indices) > 0:
                    return False
                else:
                    return True
            elif (
                len(np.where(d < 0)[0]) > 0
                and len(np.where(m[:, si_indices] == 1)[0]) > 0
            ):
                if not self.is_protection():
                    return True
                else:
                    return False
            elif len(np.where(d > 0)[0]) > 0 and len(np.where(d < 0)[0]) == 0:
                if len(indices) == 0 and len(s_indices) > 0 and len(o_indices) > 0:
                    return True
            elif (
                len(np.where(d > 0)[0]) == 0
                and len(np.where(d < 0)[0]) > 0
                and self.num_reactants == 2
                and len(x_indices) > 0
            ):
                for x_idx in x_indices:
                    if len(np.where(m[:, x_idx] < 0)[0]) == len(
                        np.where(m[:, x_idx] > 0)[0]
                    ):
                        return True
            else:
                return False
        else:
            return False
        return False

    def is_aromatic_heterocycle(self) -> bool:
        """Assesses whether the reaction involves the formation or modification of an aromatic heterocycle.

        Returns:
            bool: True if the reaction pertains to aromatic heterocycle changes, otherwise False.
        """
        aromatic_changes = np.array(
            list(np.where(self.sanitized_transformation_matrix == 1.5)[0])
            + list(np.where(self.sanitized_transformation_matrix == 0.5)[0])
        )
        if len(aromatic_changes) == 0:
            return False
        else:
            more_rings = self.ring_change > 0
            if more_rings and self.nos_reaction_center:
                a = np.array(
                    [
                        self.atom_mapping_index[atom_map]
                        for atom_map in self.transformation_mapping
                    ]
                )
                if len(a) == 0:
                    return False
                aromatic_bonds = len(np.where(self.be_matrix_products[:, a] == 1.5)[0])
                if aromatic_bonds > 0:
                    return True
                else:
                    return False
            elif self.ring_change == 0 and self.nos_reaction_center:
                a = np.array(
                    [
                        self.atom_mapping_index[atom_map]
                        for atom_map in self.transformation_mapping
                    ]
                )
                if len(a) == 0:
                    return False
                aromatic_bonds_products = len(
                    np.where(self.be_matrix_products[:, a] == 1.5)[0]
                )
                aromatic_bonds_reactants = len(
                    np.where(self.be_matrix_reactants[:, a] == 1.5)[0]
                )
                if aromatic_bonds_reactants == 0 and aromatic_bonds_products > 0:
                    return True
                else:
                    return False
        return False

    def is_reduction(self) -> bool:
        """Determines if the reaction is a reduction process based on the change in oxidation states and functional group transformation.

        Returns:
            bool: True if the reaction can be classified as a reduction, otherwise False.
        """
        if self.num_reactants == 1:
            m = self.transformation_matrix
            only_nonpositive = (
                len(np.where(self.sanitized_transformation_matrix > 0)[0]) == 0
            )
            a = np.array(
                [
                    self.atom_mapping_index[atom_map]
                    for atom_map in self.transformation_mapping
                ]
            )
            if len(a) == 0:
                return False
            involved_types = self.atoms_diagonal[a]
            o_indices = np.where(involved_types == 8)[0]
            s_indices = np.where(involved_types == 16)[0]
            if not only_nonpositive:
                n_indices = np.where(involved_types == 7)[0]
                if (
                    len(np.where(m == -0.5)[0]) > 0
                    and len(np.where(m == 0.5)[0]) > 0
                    and len(np.where(m == -1)[0]) > 0
                ):
                    return True
                elif len(s_indices) > 0 and len(o_indices) > 0:
                    if (
                        len(np.where(m.diagonal()[s_indices] == 4)[0]) > 0
                        and len(np.where(m.diagonal()[o_indices] == -4)[0]) > 0
                    ):
                        return True
                    else:
                        return False
                elif len(n_indices) > 0 and len(o_indices) > 0:  # Oxime reduction
                    if (
                        len(np.where(m.diagonal()[n_indices] == -2)[0]) > 0
                        and len(np.where(m.diagonal()[o_indices] == 0)[0]) > 0
                        and len(np.where(m[n_indices] >= 0)[0]) == 0
                    ):
                        return True
                    elif (
                        0
                        < len(np.where(m.diagonal()[o_indices] == -5)[0])
                        == len(np.where(m > 0)[0])
                        and len(np.where(m.diagonal()[n_indices] == 1)[0]) > 0
                    ):
                        return True
                    else:
                        return False
                elif len(n_indices) == 0:
                    return False
                elif (
                    len(
                        np.where(self.transformation_matrix.diagonal()[n_indices] == 1)[
                            0
                        ]
                    )
                    == 0
                ):
                    return False
                if len(np.where(self.sanitized_transformation_matrix > 0)[0]) > len(
                    np.where(
                        self.sanitized_transformation_matrix.diagonal()[n_indices] == 1
                    )[0]
                ):
                    return False
            c_indices = np.where(involved_types == 6)[0]
            if len(c_indices) == 0:
                return True
            carbons = self.mapping_diagonal[a[c_indices]]
            oxygens = self.mapping_diagonal[a[o_indices]]
            sulfurs = self.mapping_diagonal[a[s_indices]]
            all_reactant_mappings = np.array(list(self.reactant_map_dict.keys()))
            all_product_mappings = np.array(list(self.product_map_dict.keys()))
            carbon_balance_reactants = len(
                all_reactant_mappings[np.in1d(all_reactant_mappings, carbons)]
            ) == len(carbons)
            carbon_balance_products = len(
                all_product_mappings[np.in1d(all_product_mappings, carbons)]
            ) == len(carbons)
            if carbon_balance_products and carbon_balance_reactants:
                return True
            elif len(
                all_reactant_mappings[np.in1d(all_reactant_mappings, oxygens)]
            ) > len(all_product_mappings[np.in1d(all_product_mappings, oxygens)]):
                return True
            elif len(
                all_reactant_mappings[np.in1d(all_reactant_mappings, sulfurs)]
            ) > len(all_product_mappings[np.in1d(all_product_mappings, sulfurs)]):
                return True
            else:
                return False
        else:
            return False

    def is_oxidation(self) -> bool:
        """Checks if the reaction is an oxidation by examining changes in oxidation states and the involvement of key functional groups.

        Returns:
            bool: True if the reaction involves oxidation, otherwise False.
        """
        m = self.transformation_matrix
        sm = self.sanitized_transformation_matrix
        if self.num_reactants == 1:
            negative = len(np.where(m < 0)[0]) > 0
            a = np.array(
                [
                    self.atom_mapping_index[atom_map]
                    for atom_map in self.transformation_mapping
                ]
            )
            if len(a) == 0:
                return False
            involved_types = self.atoms_diagonal[a]
            if negative:
                n_indices = np.where(involved_types == 7)[0]
                s_indices = np.where(involved_types == 16)[0]
                o_indices = np.where(involved_types == 8)[0]
                c_indices = np.where(involved_types == 6)[0]
                if (
                    len(np.where(m == -0.5)[0]) > 0
                    and len(np.where(m == 0.5)[0]) > 0
                    and len(np.where(m == 1)[0]) > 0
                ):
                    return True
                elif len(s_indices) > 0 and len(o_indices) > 0:
                    if (
                        len(np.where(m.diagonal()[s_indices] == -4)[0]) > 0
                        and len(np.where(m.diagonal()[o_indices] == 4)[0]) > 0
                    ):
                        return True
                    elif (
                        len(np.where(m.diagonal()[o_indices] == 4)[0]) > 0
                        and len(np.where(m.diagonal()[s_indices] == -2)[0]) > 0
                    ):
                        return True
                    else:
                        return False
                elif len(n_indices) > 0 and len(o_indices) > 0:  # Oxime oxidation
                    if (
                        len(np.where(m.diagonal()[n_indices] == 2)[0]) > 0
                        and len(np.where(m.diagonal()[o_indices] == 0)[0]) > 0
                        and len(np.where(m[n_indices] <= 0)[0]) == 0
                    ):
                        return True
                    else:
                        return False
                elif len(c_indices) > 0 and len(o_indices) > 0:
                    try:
                        if (
                            len(np.where(m[c_indices, o_indices] < 0)[0]) == 0
                            and len(np.where(m.diagonal()[o_indices] == 4)[0]) > 0
                        ):
                            return True
                    except IndexError:
                        return False
                elif len(n_indices) == 0:
                    return False
                elif len(np.where(m.diagonal()[n_indices] == -1)[0]) == 0:
                    return False
                if len(np.where(m < 0)[0]) > len(
                    np.where(m.diagonal()[n_indices] == 1)[0]
                ):
                    return False
            d = m.diagonal()
            zero_diagonal = ~d.any()
            if not zero_diagonal:
                indices = np.where(self.r_matrix.diagonal() != 0)[0]
                ox_types = self.atoms_diagonal[indices]
                for atom_type in ox_types:
                    if atom_type != 8:
                        return False
            c_indices = np.where(involved_types == 6)[0]
            if len(c_indices) == 0:
                return True
            carbons = self.mapping_diagonal[a[c_indices]]
            all_reactant_mappings = np.array(list(self.reactant_map_dict.keys()))
            all_product_mappings = np.array(list(self.product_map_dict.keys()))
            carbon_balance_reactants = len(
                all_reactant_mappings[np.in1d(all_reactant_mappings, carbons)]
            ) == len(carbons)
            carbon_balance_products = len(
                all_product_mappings[np.in1d(all_product_mappings, carbons)]
            ) == len(carbons)
            if carbon_balance_reactants and carbon_balance_products:
                return True
            else:
                return False
        else:
            negative = len(np.where(m < 0)[0]) > 0
            a = np.array(
                [
                    self.atom_mapping_index[atom_map]
                    for atom_map in self.transformation_mapping
                ]
            )
            if len(a) == 0:
                return False
            involved_types = self.atoms_diagonal[a]
            if negative:
                s_indices = np.where(involved_types == 16)[0]
                o_indices = np.where(involved_types == 8)[0]
                sulfurs = a[s_indices]
                for sulfur in sulfurs:
                    if len(
                        np.where(self.be_matrix_products[sulfur, a[o_indices]] == 2)[0]
                    ) > len(
                        np.where(self.be_matrix_reactants[sulfur, a[o_indices]] == 2)[0]
                    ):
                        return True
            else:
                return False
            return False

    def is_acylation(self) -> bool:
        """Evaluates whether the reaction involves acylation, specifically focusing on the transformation around carbonyl groups.

        Returns:
            bool: True if the reaction involves acylation, otherwise False.
        """
        if self.num_reactants == 1:
            return False
        else:
            carbonyl_c = -1
            acyl = False
            for idx in self.reaction_center_idx:
                mp = Chem.MolFromSmiles(self.products)
                atom = mp.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() == 6:
                    for nb in atom.GetNeighbors():
                        if nb.GetAtomicNum() == 8:
                            b = mp.GetBondBetweenAtoms(idx, nb.GetIdx())
                            bo = b.GetBondTypeAsDouble()
                            if bo == 2:
                                atom_map = atom.GetAtomMapNum()
                                if atom_map in self.atom_mapping_index:
                                    carbonyl_c = self.atom_mapping_index[atom_map]
                                    acyl = True
                                break
                        else:
                            continue
                elif atom.GetAtomicNum() == 8:
                    for nb in atom.GetNeighbors():
                        if nb.GetAtomicNum() == 6:
                            for nbb in nb.GetNeighbors():
                                if nbb.GetIdx() == idx:
                                    continue
                                else:
                                    b = mp.GetBondBetweenAtoms(
                                        nbb.GetIdx(), nb.GetIdx()
                                    )
                                    bo = b.GetBondTypeAsDouble()
                                    if bo == 2:
                                        atom_map = nb.GetAtomMapNum()
                                        if atom_map in self.atom_mapping_index:
                                            carbonyl_c = self.atom_mapping_index[
                                                atom_map
                                            ]
                                            acyl = True
                                        break
                        else:
                            continue
                if acyl:
                    carbamate = mp.GetSubstructMatches(
                        Chem.MolFromSmarts(
                            "[O;H0;D2;+0]-[C;H0;D3;+0](=[O;H0;D1;+0])-[NX3;+0]"
                        )
                    )
                    rcid = np.array(self.reaction_center_idx)
                    for tup in carbamate:
                        cid = np.array(list(tup))
                        carbamate_in_rc = rcid[np.in1d(rcid, cid)]
                        if len(carbamate_in_rc) > 0:
                            return False
                        else:
                            continue
                    ester_reactant = self.mol_reactant.GetSubstructMatches(
                        Chem.MolFromSmarts(
                            "[OX2;+0]-[C;H0;D3;+0](=[O;H0;D1;+0])-[#6;!H3]"
                        )
                    )
                    for tup in ester_reactant:
                        cid = np.array(list(tup))
                        carboxyl_in_rcr = rcid[np.in1d(rcid, cid)]
                        ester_products = mp.GetSubstructMatches(
                            Chem.MolFromSmarts(
                                "[OX2;+0]-[C;H0;D3;+0](=[O;H0;D1;+0])-[#6;!H3]"
                            )
                        )
                        for tup in ester_products:
                            cid = np.array(list(tup))
                            carboxyl_in_rcp = rcid[np.in1d(rcid, cid)]
                            if len(carboxyl_in_rcr) > 0 and len(carboxyl_in_rcp) > 0:
                                if (
                                    len(
                                        np.intersect1d(carboxyl_in_rcr, carboxyl_in_rcp)
                                    )
                                    > 0
                                ):
                                    return False

            a = np.array(
                [
                    self.atom_mapping_index[atom_map]
                    for atom_map in self.transformation_mapping
                ]
            )
            if len(a) == 0:
                return False
            involved_types = self.atoms_diagonal[a]
            s_indices = np.where(involved_types == 16)[0]
            o_indices = np.where(self.atoms_diagonal == 8)[0]
            x_indices = np.array(
                list(np.where(self.atoms_diagonal == 9)[0])
                + list(np.where(self.atoms_diagonal == 17)[0])
                + list(np.where(self.atoms_diagonal == 35)[0])
                + list(np.where(self.atoms_diagonal == 53)[0])
            )
            if len(s_indices) > 0 and not self.is_oxidation():
                sulfurs = a[s_indices]
                for sulfur in sulfurs:
                    sulfone_r = np.where(
                        self.be_matrix_reactants[sulfur, o_indices] == 0
                    )[0]
                    sulfone_p = np.where(
                        self.be_matrix_products[sulfur, o_indices] == 2
                    )[0]
                    if len(sulfone_r) > 0 and len(sulfone_p) > 0:
                        if len(np.intersect1d(sulfone_r, sulfone_p)) > 0:
                            return True
                        elif len(x_indices) > 0:
                            if (
                                len(
                                    np.where(
                                        self.be_matrix_reactants[sulfur, x_indices] == 1
                                    )[0]
                                )
                                > 0
                                and len(
                                    np.where(
                                        self.be_matrix_products[sulfur, x_indices] == 1
                                    )[0]
                                )
                                == 0
                            ):
                                return True

            if acyl and carbonyl_c != -1:
                m = self.r_matrix
                bonds_broken = np.where(m[:, carbonyl_c] == -1)[0]
                bonds_formed = np.where(m[:, carbonyl_c] == 1)[0]
                bond_formed_with = []
                for bond in bonds_formed:
                    bond_formed_with.append(self.atoms_diagonal[bond])
                if (
                    7 in bond_formed_with
                    or 8 in bond_formed_with
                    or 16 in bond_formed_with
                ):
                    return True
                elif len(bonds_broken) == len(bonds_formed) and len(bonds_formed) == 1:
                    if self.be_matrix_products[carbonyl_c, bonds_broken] == 0:
                        return False
                    else:
                        return True
                else:
                    return False

            return acyl

    def is_heteroatom_alkylation(self) -> bool:
        """Determines if the reaction involves alkylation of heteroatoms (N, O, S).

        Returns:
            bool: True if the reaction is a heteroatom alkylation, otherwise False.
        """
        if not self.nos_reaction_center:
            return False
        elif self.num_reactants == 1 and self.ring_change <= 0:
            return False
        else:
            a = np.array(
                [
                    self.atom_mapping_index[atom_map]
                    for atom_map in self.transformation_mapping
                ]
            )
            if len(a) == 0:
                return False
            involved_types = self.atoms_diagonal[a]
            nos_indices = (
                list(np.where(involved_types == 7)[0])
                + list(np.where(involved_types == 8)[0])
                + list(np.where(involved_types == 16)[0])
            )
            if len(nos_indices) == 0:
                return False
            heteroatoms = a[nos_indices]
            c_indices = np.where(involved_types == 6)[0]
            o_indices = np.where(self.atoms_diagonal == 8)[0]
            if len(c_indices) == 0:
                return False
            carbons = a[c_indices]
            for heteroatom in heteroatoms:
                maps_1 = np.array(
                    list(np.where(self.r_matrix[:, heteroatom] == 1)[0])
                    + list(np.where(self.r_matrix[:, heteroatom] == 2)[0])
                )
                if len(maps_1) == 0:
                    continue
                else:
                    no_bonds = np.where(
                        self.be_matrix_reactants[maps_1, heteroatom] == 0
                    )[0]
                    carbon_bonds = maps_1[np.in1d(maps_1, carbons)]
                    # rcid = np.array(self.reaction_center_idx) <--- seems to be unused
                    carbonyls_r: list[int] = []
                    carbonyls_p: list[int] = []
                    if len(o_indices) != 0:
                        for carbon in carbons:
                            carbonyls_r += list(
                                np.where(
                                    self.be_matrix_reactants[carbon, o_indices] == 2
                                )[0]
                            )
                            carbonyls_p += list(
                                np.where(
                                    self.be_matrix_products[carbon, o_indices] == 2
                                )[0]
                            )
                    if len(carbon_bonds) > 0 and len(no_bonds) > 0:
                        if len(
                            np.intersect1d(maps_1[no_bonds], carbon_bonds)
                        ) > 0 and not (len(carbonyls_r) > 0 and len(carbonyls_p) > 0):
                            return True
                        else:
                            return False
                    else:
                        return False
            return False

    def is_cc_coupling(self) -> bool:
        """Checks if the reaction is a carbon-carbon coupling process.

        Returns:
            bool: True if the reaction involves carbon-carbon coupling, otherwise False.
        """
        if self.num_reactants == 1 and self.ring_change <= 0:
            return False
        else:
            a = np.array(
                [
                    self.atom_mapping_index[atom_map]
                    for atom_map in self.transformation_mapping
                ]
            )
            if len(a) == 0:
                return False
            involved_types = self.atoms_diagonal[a]
            c_indices = np.where(involved_types == 6)[0]
            if len(c_indices) == 0:
                return False
            carbons = a[c_indices]
            for carbon in carbons:
                maps_1 = np.array(
                    list(np.where(self.r_matrix[:, carbon] == 1)[0])
                    + list(np.where(self.r_matrix[:, carbon] == 2)[0])
                )
                if len(maps_1) == 0:
                    continue
                else:
                    no_bonds = np.where(self.be_matrix_reactants[maps_1, carbon] == 0)[
                        0
                    ]
                    if len(maps_1[np.in1d(maps_1, carbons)]) > 0 and len(no_bonds) > 0:
                        return True
                    else:
                        continue
            return False

    def is_fga(self) -> bool:
        """Determines if the reaction involves the addition of functional groups to the existing molecular framework.

        Returns:
            bool: True if the reaction is classified as functional group addition, otherwise False.
        """
        m = self.transformation_matrix
        sm = self.sanitized_transformation_matrix
        a = np.array(
            [
                self.atom_mapping_index[atom_map]
                for atom_map in self.transformation_mapping
            ]
        )
        if len(a) == 0:
            return False
        involved_types = self.atoms_diagonal[a]
        d = m.diagonal()
        halogens = np.array([9, 17, 35, 53])
        x_indices = np.where(np.in1d(self.reaction_center_atoms, halogens))[0]
        if len(np.where(d < 0)[0]) > 0:
            # Atoms that are not found in products
            count_x_additions = 0
            for x in x_indices:
                if sum(m[:, x]) == 0:
                    count_x_additions += 1

            if count_x_additions > 0 and 2 * count_x_additions == len(
                np.where(m > 0)[0]
            ):
                return True

            si_indices = np.where(involved_types == 14)[0]
            for si in si_indices:
                if (
                    len(np.where(m[:, si] == 1)[0]) > 0
                    and len(np.where(sm.diagonal() < 0)[0]) == 0
                ):
                    if not self.is_protection():
                        return True
                    else:
                        return False

            n_indices = np.where(involved_types == 7)[0]
            o_indices = np.where(involved_types == 6)[0]
            c_indices = np.where(involved_types == 8)[0]
            nitro = Chem.MolFromSmarts("[NX3;+1](=[O;H0;D1;+0])[O;H0;D1;-1]")
            if (
                len(self.mol_reactant.GetSubstructMatches(nitro)) > 0
                and len(self.mol_product.GetSubstructMatches(nitro)) > 0
                and len(n_indices) > 0
                and len(c_indices) > 0
                and self.num_reactants == 2
            ):
                for n in n_indices:
                    if len(np.where(m[:, n] == 1)[0]) > 0:
                        return True
            return False
        elif not d.any():
            metals = np.array([3, 5, 11, 12, 29, 30, 34, 47, 50])
            metal_indices = np.where(np.in1d(self.reaction_center_atoms, metals))[0]
            if len(metal_indices) > 0:
                for metal_idx in metal_indices:
                    if len(np.where(m[:, metal_idx] < 0)[0]) == 0:
                        return True
                    else:
                        return False
            elif len(x_indices) > 0 and len(np.where(m < 0)[0]) == 0:
                return True
            else:
                return False
        else:
            c_indices = np.where(involved_types == 6)[0]
            if len(c_indices) == 0:
                return False
            carbons = a[c_indices]
            for carbon in carbons:
                maps_1 = np.array(
                    list(np.where(self.r_matrix[:, carbon] == 1)[0])
                    + list(np.where(self.r_matrix[:, carbon] == 2)[0])
                )
                if len(maps_1) == 0:
                    continue
                else:
                    no_bonds = np.where(self.be_matrix_reactants[maps_1, carbon] == 0)[
                        0
                    ]
                    if len(maps_1[np.in1d(maps_1, carbons)]) > 0 and len(no_bonds) > 0:
                        continue
                    else:
                        return True

            return False
        return False

    def is_deprotection(self) -> bool:
        """Evaluates whether the reaction is a deprotection, which involves the removal of protective groups from functional sites.

        Returns:
            bool: True if the reaction is a deprotection process, otherwise False.
        """
        if self.num_reactants == 1:
            only_nonpositive = (
                len(np.where(self.sanitized_transformation_matrix > 0)[0]) == 0
            )
            if not only_nonpositive:
                return False
            a = np.array(
                [
                    self.atom_mapping_index[atom_map]
                    for atom_map in self.transformation_mapping
                ]
            )
            if len(a) == 0:
                return False
            involved_types = self.atoms_diagonal[a]
            nos_indices = (
                list(np.where(involved_types == 7)[0])
                + list(np.where(involved_types == 8)[0])
                + list(np.where(involved_types == 16)[0])
            )
            c_si_indices = list(np.where(involved_types == 6)[0]) + list(
                np.where(involved_types == 14)[0]
            )
            if len(nos_indices) == 0:
                alkyne_indices = []
                for carbon in c_si_indices:
                    c_idx = a[carbon]
                    if len(np.where(self.be_matrix_reactants[:, c_idx] == 3)[0]) > 0:
                        alkyne_indices.append(carbon)
                if len(alkyne_indices) == 0:
                    return False
                else:
                    nos_indices += alkyne_indices
            heteroatoms = a[nos_indices]
            if len(c_si_indices) == 0:
                return False
            carbons = a[c_si_indices]
            for heteroatom in heteroatoms:
                maps_1 = np.where(self.r_matrix[:, heteroatom] == -1)[0]
                if len(maps_1) == 0:
                    continue
                else:
                    no_bonds = len(
                        np.where(self.be_matrix_products[maps_1, heteroatom] == 0)[0]
                    )
                    if len(maps_1[np.in1d(maps_1, carbons)]) > 0 and no_bonds > 0:
                        in_products = len(
                            np.where(
                                self.be_matrix_products[
                                    maps_1[np.in1d(maps_1, carbons)]
                                ]
                                != 0
                            )[0]
                        )
                        if in_products > 0:
                            continue
                        else:
                            return True
                    else:
                        continue
            return False

        else:
            return False

    def is_protection(self) -> bool:
        """Determines if the reaction is a protection, which involves adding protective groups to functional sites.

        Returns:
            bool: True if the reaction is classified as a protection process, otherwise False.
        """
        if self.num_reactants == 1:
            return False
        else:
            a = np.array(
                [
                    self.atom_mapping_index[atom_map]
                    for atom_map in self.transformation_mapping
                ]
            )
            if len(a) == 0:
                return False
            involved_types = self.atoms_diagonal[a]
            nos_indices = (
                list(np.where(involved_types == 7)[0])
                + list(np.where(involved_types == 8)[0])
                + list(np.where(involved_types == 16)[0])
            )
            c_indices = list(np.where(involved_types == 6)[0])
            si_indices = list(np.where(involved_types == 14)[0])
            c_si_indices = c_indices + si_indices
            alkyne_indices = []
            if len(nos_indices) == 0:
                for carbon in c_indices:
                    c_idx = a[carbon]
                    if len(np.where(self.be_matrix_reactants[:, c_idx] == 3)[0]) > 0:
                        alkyne_indices.append(carbon)
                if len(alkyne_indices) == 0:
                    return False
                else:
                    nos_indices += alkyne_indices
            heteroatoms = a[nos_indices]
            if len(c_si_indices) == 0:
                return False
            carbons = a[c_si_indices]
            if len(si_indices) > 0:
                silicons = a[si_indices]
            else:
                silicons = []
            for heteroatom in heteroatoms:
                maps_1 = np.where(self.r_matrix[:, heteroatom] == 1)[0]
                if len(maps_1) == 0:
                    continue
                elif heteroatom in carbons:  # Indicates an alkyne
                    if len(silicons) == 0:
                        continue
                    elif len(maps_1[np.in1d(maps_1, silicons)]) > 0:
                        return True
                    else:
                        continue
                else:
                    no_bonds = len(
                        np.where(self.be_matrix_reactants[maps_1, heteroatom] == 0)[0]
                    )
                    if len(maps_1[np.in1d(maps_1, carbons)]) > 0 and no_bonds > 0:
                        return True
                    else:
                        continue
            return False

    def classify_reaction(self) -> str:
        """Classifies the reaction based on its chemical characteristics and transformation patterns.

        Returns:
            str: The classification of the reaction, such as 'Reduction', 'Oxidation', etc.
        """
        if self.is_aromatic_heterocycle():
            return "Aromatic Heterocycle Formation"
        elif self.is_acylation():
            return "Acylation"
        elif self.is_fgi():
            return "Functional Group Interconversion"
        elif self.is_reduction():
            return "Reduction"
        elif self.is_oxidation():
            return "Oxidation"
        elif self.is_fga():
            return "Functional Group Addition"
        elif self.is_heteroatom_alkylation():
            return "Heteroatom Alkylation and Arylation"
        elif self.is_cc_coupling():
            return "C-C Coupling"
        elif self.is_deprotection():
            return "Deprotection"
        elif self.is_protection():
            return "Protection"
        else:
            return "Miscellaneous"

    def name_reaction(self, smirks_db: pd.DataFrame) -> str:
        """Determines the name of the reaction from a database based on SMIRKS transformations.

        Args:
            smirks_db (pd.DataFrame): DataFrame containing SMIRKS patterns and corresponding reaction names.

        Returns:
            str: The name of the reaction, or 'OtherReaction' if no specific name can be determined.
        """
        reactants_smiles, products_smiles = self.sanitized_reaction.split(">>")
        reactants = reactants_smiles.split(".")
        products = products_smiles.split(".")

        if (
            len(reactants) > 4 or len(products) > 4
        ):  # There are no templates for reactions with more than four reactants.
            return "OtherReaction"

        new_products = []  # Try to canonicalize SMILES

        for product in products:
            try:
                new_products.append(
                    Chem.MolToSmiles(Chem.MolFromSmiles(product), isomericSmiles=False)
                )
            except:
                new_products.append(product)

        num_reactants = len(reactants)
        # num_products = len(products)

        rxn_name = ""
        selected_rxns = smirks_db[smirks_db["nreact"] == num_reactants]
        react_tuple = tuple(Chem.MolFromSmiles(reactant) for reactant in reactants)

        if num_reactants == 1:
            all_tuples = [react_tuple]
        else:
            all_tuples = list(
                itertools.permutations(react_tuple)
            )  # RDKit does not permute reactants by itself

        # TODO: Further refine reactions by superclass

        for i in selected_rxns.index:  # Iterate over all reactants to find a match
            smirks = selected_rxns["smirks"][i]
            rxn = AllChem.ReactionFromSmarts(smirks)
            pred_products = []

            for tup in all_tuples:
                try:
                    pred_product = rxn.RunReactants(tup)
                except Exception:
                    continue
                pred_products += pred_product

            if len(pred_products) == 0:  # No products are found
                continue
            else:
                for prods in pred_products:
                    try:
                        prod = Chem.MolToSmiles(prods[0], isomericSmiles=False)
                    except Exception:
                        continue

                    if (
                        prod in new_products
                    ):  # Predicted product is in the real reaction
                        rxn_name = selected_rxns["name"][i].strip("{}")
                        return rxn_name
                    else:
                        continue

        if rxn_name == "":
            rxn_name = "OtherReaction"

        return rxn_name


if __name__ == "__main__":
    rxn_smiles_with_atom_mapping = "[CH3:9][CH:8]([CH3:10])[c:7]1[cH:11][cH:12][cH:13][cH:14][cH:15]1.[CH2:3]([CH2:4][C:5](=[O:6])Cl)[CH2:2][Cl:1].[Al+3].[Cl-].[Cl-].[Cl-].C(Cl)Cl>>[CH3:9][CH:8]([CH3:10])[c:7]1[cH:11][cH:12][c:13]([cH:14][cH:15]1)[C:5](=[O:6])[CH2:4][CH2:3][CH2:2][Cl:1]"

    ReactionClassifier(rxn_smiles_with_atom_mapping, keep_mapping=True)
