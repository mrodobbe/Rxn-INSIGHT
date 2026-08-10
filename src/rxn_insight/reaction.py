"""Reaction module"""

from __future__ import annotations

import json
import hashlib
import warnings
from typing import TYPE_CHECKING, Optional, Union, Dict, Any

import numpy as np
import pandas as pd
from rdkit import Chem
from tqdm import tqdm

from rxn_insight.classification import ReactionClassifier
from rxn_insight.naming import name_reaction
from rxn_insight.utils import (
    curate_smirks,
    draw_chemical_reaction,
    get_catalyst_ranking,
    get_fp,
    get_reaction_template,
    get_reagent_ranking,
    get_scaffold,
    get_similarity,
    get_solvent_ranking,
    remove_atom_mapping,
)
if TYPE_CHECKING:  # imported lazily at runtime to keep `import rxn_insight` fast
    from rxnmapper import RXNMapper

import logging
logger = logging.getLogger(__name__)


def _element_symbol(atomic_number) -> str:
    """Return the element symbol for an atomic number (RDKit periodic table)."""
    return Chem.GetPeriodicTable().GetElementSymbol(int(atomic_number))


class Reaction:

    """Handles operations related to chemical reactions.

    This class facilitates various operations on chemical reactions, such as
    parsing reaction strings, identifying components like solvents and reagents,
    classifying reactions, and analyzing ring structures.

    Attributes:
        reaction (str): The SMILES representation of the reaction.
        solvent (str): Solvents used in the reaction.
        reagent (str): Reagents used in the reaction.
        catalyst (str): Catalysts used in the reaction.
        reference (str): Reference or note associated with the reaction.
        smirks_db (pd.DataFrame): Database of SMIRKS transformations.
        fg_db (pd.DataFrame): Functional group data.
        classifier (ReactionClassifier): Reaction classification object.
        reactants (str): SMILES string of the reactants.
        products (str): SMILES string of the products.
        mapped_reaction (str): Reaction with atom mappings included.
        reaction_class (str): Class of the reaction.
        template (str): Reaction template derived from the classifier via RDChiral.
        reaction_info (dict): Additional information about the reaction.
        tag (str): Optional tag for the reaction.
        name (str): Optional name of the reaction.
        byproducts (tuple): Tuple of byproducts in the reaction.
        scaffold (str): Molecular scaffold of the reaction.
        neighbors (Any): Placeholder for reaction neighborhood information.
        suggested_solvent (str): Suggested solvent for the reaction.
        suggested_catalyst (str): Suggested catalyst for the reaction.
        suggested_reagent (str): Suggested reagent for the reaction.

    Example:
        >>> import rxn_insight as ri
        >>> rxn = ri.Reaction("OB(O)c1ccccc1.Brc1ccccc1>>c1ccc(-c2ccccc2)cc1")
        >>> ri = rxn.get_reaction_info()
        >>> print(ri)
        {'REACTION': 'Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1',
        'MAPPED_REACTION': 'Br[c:5]1[cH:6][cH:7][cH:8][cH:9][cH:10]1.OB(O)[c:4]1[cH:3][cH:2][cH:1][cH:12][cH:11]1>>[cH:1]1[cH:2][cH:3][c:4](-[c:5]2[cH:6][cH:7][cH:8][cH:9][cH:10]2)[cH:11][cH:12]1',
        'N_REACTANTS': 2, 'N_PRODUCTS': 1, 'FG_REACTANTS': ['Aromatic halide', 'Boronic acid'], 'FG_PRODUCTS': [],
        'PARTICIPATING_RINGS_REACTANTS': ['c1ccccc1', 'c1ccccc1'], 'PARTICIPATING_RINGS_PRODUCTS': ['c1ccccc1', 'c1ccccc1'],
        'ALL_RINGS_PRODUCTS': ['c1ccccc1', 'c1ccccc1'], 'BY-PRODUCTS': ['HBr', 'B'], 'CLASS': 'C-C Coupling',
        'TAG': 'd79a78c79f0c392f0911481acf5c300cc98205269acdb93c24fb610a61c4c868', 'SOLVENT': [''], 'REAGENT': [''],
        'CATALYST': [''], 'REF': '', 'NAME': 'Suzuki coupling with boronic acids', 'SCAFFOLD': 'c1ccc(-c2ccccc2)cc1'}

    """

    def __init__(
        self,
        reaction: str,
        solvent: str = "",
        reagent: str = "",
        catalyst: str = "",
        ref: str = "",
        rxn_mapper: Optional[RXNMapper] = None,
        keep_mapping: bool = False,
        smirks: pd.DataFrame = None,
        fg: pd.DataFrame = None,
        classify: bool = True,
        search_template: bool = True,
        include_hydrogens: bool = False,
    ):

        """Initializes a Reaction object with details of the reaction.

        Args:
            reaction (str): A string representing the reaction in SMILES format.
            solvent (str, optional): Solvent(s) used in the reaction. Defaults to an empty string.
            reagent (str, optional): Reagent(s) used in the reaction. Defaults to an empty string.
            catalyst (str, optional): Catalyst(s) used in the reaction. Defaults to an empty string.
            ref (str, optional): Reference or note associated with the reaction. Defaults to an empty string.
            rxn_mapper (RXNMapper, optional): Object for reaction mapping. Defaults to None.
            keep_mapping (bool, optional): Whether to retain atom mappings in the reaction. Defaults to False.
            smirks (pd.DataFrame, optional): DataFrame of SMIRKS transformations. Defaults to None.
            fg (pd.DataFrame, optional): DataFrame of functional groups data. Defaults to None.
            search_template (bool, optional): Whether to search for reaction templates. Defaults to True.
        """

        self.reaction = ""
        self.solvent = solvent
        self.reagent = reagent
        self.catalyst = catalyst
        self.reference = ref
        self.read_reaction(reaction)
        self.classify = classify
        if ":" in self.reaction and not keep_mapping:
            self.reaction = remove_atom_mapping(
                self.reaction
            )  # Remove atom mapping for consistency
        else:
            self.reaction = self.reaction
        self.smirks_db = smirks
        self.fg_db = fg
        if self.classify:
            self.classifier = ReactionClassifier(
                reaction,
                rxn_mapper=rxn_mapper,
                keep_mapping=keep_mapping,
                search_template=search_template,
                include_hydrogens=include_hydrogens,
            )
            self.add_agents()
            self.reactants, self.products = self.classifier.sanitized_reaction.split(">>")
            self.mapped_reaction = self.classifier.sanitized_mapped_reaction
            self.template = self.classifier.template
            self.scaffold = self.get_scaffold()

        self.reaction_class = ""
        self.reaction_info: dict[str, tuple[str, ...] | str] = dict()
        self.tag = ""
        self.name = ""
        self.byproducts: tuple[str, ...] = tuple()
        self.neighbors = None
        self.suggested_solvent = ""
        self.suggested_catalyst = ""
        self.suggested_reagent = ""

    def read_reaction(
            self,
            reaction: str
    ) -> None:

        """Processes a reaction string in SMILES format.

        Args:
            reaction (str): Reaction string in SMILES format, with components separated by `>`.
        """

        reaction_elements = reaction.split(">")
        self.reaction = f"{reaction_elements[0]}>>{reaction_elements[2]}"
        reagents = reaction_elements[1].split(".")
        if len(reagents) == 1 and reagents[0] == "":
            self.reagent = ""
        else:
            solvents = self.solvent.split(".")
            catalysts = self.catalyst.split(".")
            agents = []
            for reagent in reagents:
                if reagent in solvents or reagent in catalysts:
                    continue
                else:
                    agents.append(reagent)
            self.reagent = ".".join(agents)

    def add_agents(self) -> None:
        """Adds agents identified by the classifier to the reagent list."""
        reagents = self.reagent.split(".")
        reagents += self.classifier.extra_agents
        self.reagent = ".".join(reagents)

    def get_class(self) -> str:
        """Determines and returns the class of the reaction."""
        self.reaction_class = self.classifier.classify_reaction()
        return self.reaction_class

    def get_rings_in_products(self) -> list[str]:
        """Identifies and returns ring structures in the reaction products."""
        return self.classifier.get_ring_type(self.classifier.mol_product)

    def get_rings_in_reactants(self) -> list[str]:
        """Identifies and returns ring structures in the reaction reactants."""
        return self.classifier.get_ring_type(self.classifier.mol_reactant)

    def get_rings_in_reaction_center(
        self,
    ) -> tuple[list[str], ...]:
        """Identifies and returns rings in the reaction center for reactants and products."""
        return tuple(
            [
                self.classifier.get_ring_type(
                    self.classifier.mol_reactant, self.classifier.reactant_map_dict
                ),
                self.classifier.get_ring_type(
                    self.classifier.mol_product, self.classifier.product_map_dict
                ),
            ]
        )

    def get_functional_groups(self) -> tuple[list[str], ...]:
        """Identifies and returns functional groups in reactants and products."""
        if self.fg_db is None:
            from importlib.resources import files

            path = files(f"{__package__}.data").joinpath("functional_groups.json")
            self.fg_db = pd.read_json(path, orient="records", lines=True)
        c = self.classifier
        return tuple(
            [
                c.get_functional_groups(
                    c.mol_reactant, c.reactant_map_dict, self.fg_db
                ),
                c.get_functional_groups(c.mol_product, c.product_map_dict, self.fg_db),
            ]
        )

    def get_byproducts(self) -> list[str]:
        """Calculates and returns byproducts of the reaction based on functional group analysis."""
        fg_r, fg_p = self.get_functional_groups()
        calculated_byproducts = self.classifier.balance_reaction(fg_r, fg_p)
        self.byproducts = calculated_byproducts
        return calculated_byproducts

    def get_scaffold(self) -> Optional[str]:
        """Extracts and returns the molecular scaffold of the product."""
        return get_scaffold(self.classifier.mol_product)

    def get_name(self) -> str:
        """Determines and returns the name of the reaction based on SMIRKS data."""
        if self.smirks_db is None:
            from importlib.resources import files

            path = files(f"{__package__}.data").joinpath("smirks.json")
            self.smirks_db = curate_smirks(
                pd.read_json(path, orient="records", lines=True)
            )
        if self.classify:
            self.name = name_reaction(self.classifier.sanitized_reaction, self.smirks_db)
        else:
            self.name = name_reaction(self.reaction, self.smirks_db)

        return self.name

    def get_reaction_info(self) -> dict[str, list[str] | str]:

        """This function compiles all reaction-related information at once. Upon calling this function,
        the T-matrix of the reaction will be calculated, a class and name will be assigned, the functional groups,
        rings, and scaffold of the reaction are determined. All information is returned as a dictionary."""

        if self.fg_db is None:
            from importlib.resources import files

            path = files(f"{__package__}.data").joinpath("functional_groups.json")
            self.fg_db = pd.read_json(path, orient="records", lines=True)

        if self.classify:
            info_dict = self.classifier.get_reaction_center_info(self.fg_db)
            self.tag = info_dict["TAG"]
            self.reaction_class = info_dict["CLASS"]
        else:
            info_dict = {"REACTION": self.reaction, "NAME": self.get_name()}

        try:
            info_dict["SOLVENT"] = self.solvent.split(".")
        except AttributeError:
            info_dict["SOLVENT"] = []
        try:
            info_dict["REAGENT"] = self.reagent.split(".")
        except AttributeError:
            info_dict["REAGENT"] = []
        try:
            info_dict["CATALYST"] = self.catalyst.split(".")
        except AttributeError:
            info_dict["CATALYST"] = []
        try:
            info_dict["REF"] = self.reference
        except AttributeError:
            info_dict["REF"] = ""

        if self.name == "":
            info_dict["NAME"] = self.get_name()
        else:
            info_dict["NAME"] = self.name

        info_dict["SCAFFOLD"] = self.get_scaffold()

        self.name = info_dict["NAME"]
        self.scaffold = info_dict["SCAFFOLD"]

        self.reaction_info = info_dict

        return info_dict

    def get_reaction_center(self) -> Optional[str]:
        """Returns the reaction center SMILES string if available."""
        return self.classifier.template_smiles

    def find_neighbors(
        self,
        df: pd.DataFrame,
        fp: str = "MACCS",
        concatenate: bool = True,
        max_return: int = 100,
        threshold: float = 0.3,
        broaden: bool = False,
        full_search: bool = False,
    ) -> pd.DataFrame:
        """Finds and returns similar reactions in the database.

        Args:
            df: The DataFrame to search within.
            fp: The type of fingerprint to use, 'MACCS' or 'Morgan'.
            concatenate: Whether to concatenate patterns in fingerprinting.
            max_return: Maximum number of similar reactions to return.
            threshold: The similarity threshold to consider for matching.
            broaden: Whether to use a broadened search criteria based on tags.
            full_search: If true, performs an exhaustive search across the database.

        Example:
            >>> from rxn_insight.reaction import Reaction
            >>> df_uspto = pd.read_parquet("uspto_rxn_insight.gzip")  # Download: https://zenodo.org/records/10171745
            >>> rxn = Reaction("OB(O)c1ccccc1.Brc1ccccc1>>c1ccc(-c2ccccc2)cc1")
            >>> df_neighbors = rxn.find_neighbors(df_uspto)
        """
        self.get_reaction_info()
        if full_search:
            warnings.warn("Full database search is activated. This may take long.")
            df_tag = df.copy()
        elif broaden:
            tag = self.give_broad_tag()
            df_tag = df[df["TAG2"] == tag].copy()
        else:
            tag = self.tag
            df_tag = df[df["TAG"] == tag].copy()
        if len(df_tag.index) == 0:
            logger.warning("No similar reactions found...")
            return None
        fps = []
        if fp.lower() == "maccs" and concatenate:
            if "rxn_str_patt_fp" in df_tag:
                fps = [
                    np.fromiter(fp, dtype=np.int64)
                    for fp in tqdm(
                        df_tag["rxn_str_patt_fp"].tolist(),
                        desc="Loading fingerprints...",
                    )
                ]
        elif fp.lower() == "maccs" and not concatenate:
            if "rxn_dif_patt_fp" in df_tag:
                fps = [
                    np.fromiter(fp, dtype=np.int64)
                    for fp in tqdm(
                        df_tag["rxn_dif_patt_fp"].tolist(),
                        desc="Loading fingerprints...",
                    )
                ]
        elif fp.lower() == "morgan" and concatenate:
            if "rxn_str_morgan_fp" in df_tag:
                fps = [
                    np.fromiter(fp, dtype=np.int64)
                    for fp in tqdm(
                        df_tag["rxn_str_morgan_fp"].tolist(),
                        desc="Loading fingerprints...",
                    )
                ]
        elif fp.lower() == "morgan" and not concatenate:
            if "rxn_dif_morgan_fp" in df_tag:
                fps = [
                    np.fromiter(fp, dtype=np.int64)
                    for fp in tqdm(
                        df_tag["rxn_dif_morgan_fp"].tolist(),
                        desc="Loading fingerprints...",
                    )
                ]
        else:
            raise KeyError(
                f"Fingerprint choice {fp} is not supported. Select either MACCS or Morgan."
            )
        if len(fps) == 0:
            fps = [
                get_fp(r, fp, concatenate)
                for r in tqdm(
                    df_tag["REACTION"].tolist(), desc="Creating fingerprints..."
                )
            ]
        rxnfp = get_fp(self.reaction, fp, concatenate)

        sims = [
            get_similarity(rxnfp, fp)
            for fp in tqdm(fps, desc="Calculating Tanimoto similarity")
        ]
        df_tag["SIMILARITY"] = sims
        df_tag = df_tag.sort_values(by="SIMILARITY", ascending=False)
        df_tag["SOLVENT"] = df_tag["SOLVENT"].fillna("")
        df_tag["CATALYST"] = df_tag["CATALYST"].fillna("")
        df_tag["REAGENT"] = df_tag["REAGENT"].fillna("")
        max_similarity = df_tag["SIMILARITY"].max()
        df_tag = df_tag[df_tag["SIMILARITY"] > threshold].copy()
        logger.warning(
            f"Reaction found with similarity of {max_similarity:.3f}. This will be our best match."
        )
        df_return = df_tag.iloc[:max_return].copy()
        if "rxn_str_patt_fp" in df_return.keys():
            df_return = df_return.drop(columns=["rxn_str_patt_fp"])
        if "rxn_dif_patt_fp" in df_return.keys():
            df_return = df_return.drop(columns=["rxn_dif_patt_fp"])
        if "rxn_str_morgan_fp" in df_return.keys():
            df_return = df_return.drop(columns=["rxn_str_morgan_fp"])
        if "rxn_dif_morgan_fp" in df_return.keys():
            df_return = df_return.drop(columns=["rxn_dif_morgan_fp"])
        if "TAG" in df_return.keys():
            df_return = df_return.drop(columns=["TAG"])
        if "TAG2" in df_return.keys():
            df_return = df_return.drop(columns=["TAG2"])

        self.neighbors = df_return

        return df_return

    def give_broad_tag(self) -> str:
        """Generates a broadened tag for the reaction based on its characteristics."""
        rxn_info = self.reaction_info
        tag = f"{rxn_info['CLASS']} "
        try:
            fg_r = sorted(list(rxn_info["FG_REACTANTS"]))
        except AttributeError:
            fg_r = ""
        try:
            fg_p = sorted(list(rxn_info["FG_PRODUCTS"]))
        except AttributeError:
            fg_p = ""
        tag += " ".join(fg_r) + " "
        tag += " ".join(fg_p)
        tag_bytes = tag.encode("UTF-8")
        hashtag = hashlib.sha256(tag_bytes).hexdigest()
        return str(hashtag)

    def suggest_conditions(self, df: pd.DataFrame) -> dict[str, pd.DataFrame]:
        """Suggests reaction conditions based on similar reactions found.

        Args:
            df: The DataFrame containing reaction data to analyze.

        Example:
            >>> import rxn_insight as ri
            >>> df_uspto = pd.read_parquet("uspto_rxn_insight.gzip")  # Download: https://zenodo.org/records/10171745
            >>> rxn = ri.Reaction("OB(O)c1ccccc1.Brc1ccccc1>>c1ccc(-c2ccccc2)cc1")
            >>> df_conditions = rxn.suggest_conditions(df_uspto)

        """
        if self.neighbors is None or len(self.neighbors.index) == 0:
            nbs = self.find_neighbors(df, max_return=5000, threshold=0.3, broaden=True)
        else:
            nbs = self.neighbors

        solvent_rank = get_solvent_ranking(nbs)
        solvent_rank = solvent_rank.copy().sort_values(by="COUNT", ascending=False)
        catalyst_rank = get_catalyst_ranking(nbs)
        catalyst_rank = catalyst_rank.copy().sort_values(by="COUNT", ascending=False)
        reagent_rank = get_reagent_ranking(nbs)
        reagent_rank = reagent_rank.copy().sort_values(by="COUNT", ascending=False)

        conditions_dict = {
            "Solvent": solvent_rank["NAME"][solvent_rank.index[0]],
            "Catalyst": catalyst_rank["NAME"][catalyst_rank.index[0]],
            "Reagent": reagent_rank["NAME"][reagent_rank.index[0]],
        }
        self.suggested_solvent = solvent_rank
        self.suggested_catalyst = catalyst_rank
        self.suggested_reagent = reagent_rank

        return conditions_dict

    def draw(self, include_mapping: bool = False, filename: Union[str, None] = None) -> Union[str, None]:
        """
        Visualizes the chemical reaction as an SVG image.

        Args:
            include_mapping (bool): If True, displays the reaction with atom mappings.
                                    If False, displays without mappings. Default is False.
            filename (Union[str, None]): If provided, saves the SVG to this file path.
                                        Should include the .svg extension. Default is None.

        Returns:
            Union[str, None]: The SVG string if successful, None if there's an error.

        Example:
            >>> rxn = Reaction("CC(C)c1ccccc1.Cl>>CC(C)c1ccccc1Cl")
            >>> rxn.draw()  # Display in Jupyter
            >>> rxn.draw(filename="my_reaction.svg")  # Save to file
            >>> rxn.draw(include_mapping=True, filename="mapped_reaction.svg")  # Save with mappings
        """
        # Choose which reaction SMILES to use based on include_mapping parameter
        if include_mapping and hasattr(self, 'mapped_reaction'):
            reaction_to_draw = self.mapped_reaction
        else:
            reaction_to_draw = self.reaction

        # Generate the SVG
        svg_string = draw_chemical_reaction(reaction_to_draw)

        # Save to file if filename is provided
        if filename is not None:
            try:
                with open(filename, 'w', encoding='utf-8') as f:
                    f.write(svg_string)
                logger.warning(f"Reaction saved to {filename}")
            except IOError as e:
                logger.warning(f"Error saving file: {e}")
                return None

        # Try to display in IPython/Jupyter environment
        try:
            from IPython.display import SVG, display
            display(SVG(svg_string))
        except ImportError:
            if filename is None:
                logger.warning("This function requires IPython to be installed for display: pip install ipython")
                logger.warning("Alternatively, provide a filename to save the SVG.")
                return None

        return svg_string


    def explain(
            self,
            preserve_stereo: bool = True,
            include_ring_info: bool = True,
            renumber_mappings: bool = True,
    ) -> Dict[str, Any]:
        """
        Provide detailed explanation of reaction transformations and stereochemistry changes.

        This function analyzes a reaction to identify what changes occur, including
        bond changes, stereochemistry inversions, functional group transformations,
        and leaving groups.

        Args:
            preserve_stereo: Whether to analyze stereochemistry changes.
            include_ring_info: Whether to include ring membership information (default: True).
            renumber_mappings: Whether to renumber atom mappings to start from 1 (default: True).

        Returns:
            Dictionary containing detailed reaction analysis:
                - reaction_center: Set of atom map numbers involved in transformation
                - stereochemistry_changes: List of stereochemistry changes
                - bond_changes: Description of bond formations/breakages
                - leaving_groups: Identified leaving groups
                - template: Generated SMIRKS template
                - classification: Reaction class and name from Rxn-INSIGHT

        Examples:
            >>> rsmi = "C[C@]1(c2ccccc2)O[C@H]1c1ccccc1.C[Si](C)(C)Cl>>C[C@@]1(c2ccccc2)O[C@@]1(c1ccccc1)[Si](C)(C)C"
            >>> r = Reaction(rsmi)
            >>> explanation = r.explain()
            >>> print(explanation['stereochemistry_changes'])
        """
        bond_name = {1: "single", 1.5: "aromatic", 2: "double", 3: "triple"}
        reaction_info = self.get_reaction_info()
        classifier = self.classifier
        explanation: Dict[str, Any] = {
            "reaction_center": set(),
            "stereochemistry_changes": [],
            "bond_changes": [],
            "leaving_groups": [],
            "added_atoms": [],
            "rings": {},
            "functional_groups": {},
            "template": None,
            "classification": {}
        }

        # Get basic reaction center from transformation mapping
        reaction_center = set(classifier.transformation_mapping)

        # Analyze stereochemistry changes if requested
        if preserve_stereo:
            reactant_atom_info = {}
            product_atom_info = {}
            reactant_bond_info = {}  # For E/Z stereochemistry
            product_bond_info = {}   # For E/Z stereochemistry

            # Collect reactant stereochemistry
            for molr in classifier.reactant_mols:
                Chem.AssignStereochemistry(molr, cleanIt=False, force=True)
                for atom in molr.GetAtoms():
                    if atom.HasProp('molAtomMapNumber'):
                        map_num = atom.GetIntProp('molAtomMapNumber')
                        chirality = atom.GetChiralTag()
                        cip = atom.GetProp('_CIPCode') if atom.HasProp('_CIPCode') else None
                        reactant_atom_info[map_num] = {
                            'mol': molr,
                            'idx': atom.GetIdx(),
                            'chirality': chirality,
                            'cip': cip,
                            'symbol': atom.GetSymbol()
                        }

                # Collect bond stereochemistry (E/Z)
                for bond in molr.GetBonds():
                    begin_atom = bond.GetBeginAtom()
                    end_atom = bond.GetEndAtom()

                    if begin_atom.HasProp('molAtomMapNumber') and end_atom.HasProp('molAtomMapNumber'):
                        begin_map = begin_atom.GetIntProp('molAtomMapNumber')
                        end_map = end_atom.GetIntProp('molAtomMapNumber')
                        bond_key = tuple(sorted([begin_map, end_map]))

                        stereo = bond.GetStereo()
                        bond_type = bond.GetBondTypeAsDouble()

                        reactant_bond_info[bond_key] = {
                            'bond_type': bond_type,
                            'stereo': stereo,
                            'atoms': (begin_map, end_map)
                        }

            # Collect product stereochemistry
            for molp in classifier.product_mols:
                Chem.AssignStereochemistry(molp, cleanIt=False, force=True)
                for atom in molp.GetAtoms():
                    if atom.HasProp('molAtomMapNumber'):
                        map_num = atom.GetIntProp('molAtomMapNumber')
                        chirality = atom.GetChiralTag()
                        cip = atom.GetProp('_CIPCode') if atom.HasProp('_CIPCode') else None
                        product_atom_info[map_num] = {
                            'mol': molp,
                            'idx': atom.GetIdx(),
                            'chirality': chirality,
                            'cip': cip,
                            'symbol': atom.GetSymbol()
                        }

                # Collect bond stereochemistry (E/Z)
                for bond in molp.GetBonds():
                    begin_atom = bond.GetBeginAtom()
                    end_atom = bond.GetEndAtom()

                    if begin_atom.HasProp('molAtomMapNumber') and end_atom.HasProp('molAtomMapNumber'):
                        begin_map = begin_atom.GetIntProp('molAtomMapNumber')
                        end_map = end_atom.GetIntProp('molAtomMapNumber')
                        bond_key = tuple(sorted([begin_map, end_map]))

                        stereo = bond.GetStereo()
                        bond_type = bond.GetBondTypeAsDouble()

                        product_bond_info[bond_key] = {
                            'bond_type': bond_type,
                            'stereo': stereo,
                            'atoms': (begin_map, end_map)
                        }

            # Identify stereochemistry changes
            for map_num in reactant_atom_info:
                if map_num in product_atom_info:
                    r_info = reactant_atom_info[map_num]
                    p_info = product_atom_info[map_num]

                    if r_info['chirality'] != p_info['chirality'] or r_info['cip'] != p_info['cip']:
                        reaction_center.add(map_num)

                        change_desc = {
                            'atom': f"{r_info['symbol']}:{map_num}",
                            'chirality_change': None,
                            'cip_change': None
                        }

                        if r_info['chirality'] != p_info['chirality']:
                            change_desc['chirality_change'] = f"{r_info['chirality']} -> {p_info['chirality']}"

                        if r_info['cip'] != p_info['cip']:
                            change_desc['cip_change'] = f"{r_info['cip']} -> {p_info['cip']}"

                        explanation['stereochemistry_changes'].append(change_desc)

            # Identify E/Z stereochemistry changes
            all_bond_keys = set(reactant_bond_info.keys()) | set(product_bond_info.keys())

            for bond_key in all_bond_keys:
                r_bond = reactant_bond_info.get(bond_key, None)
                p_bond = product_bond_info.get(bond_key, None)

                # Case 1: Triple bond -> Stereo double bond (creation of E/Z)
                if r_bond and p_bond:
                    if r_bond['bond_type'] == 3.0 and p_bond['bond_type'] == 2.0:
                        if p_bond['stereo'] in (Chem.BondStereo.STEREOZ, Chem.BondStereo.STEREOE):
                            atom1, atom2 = bond_key
                            sym1 = reactant_atom_info.get(atom1, {}).get('symbol', 'X')
                            sym2 = reactant_atom_info.get(atom2, {}).get('symbol', 'X')

                            config = "Z (cis)" if p_bond['stereo'] == Chem.BondStereo.STEREOZ else "E (trans)"
                            change_desc = {
                                'bond': f"{sym1}:{atom1}={sym2}:{atom2}",
                                'type': 'E/Z',
                                'change': f"Triple bond reduced to {config} double bond"
                            }
                            explanation['stereochemistry_changes'].append(change_desc)
                            reaction_center.add(atom1)
                            reaction_center.add(atom2)

                    # Case 2: Single bond -> Stereo double bond (creation of E/Z)
                    elif r_bond['bond_type'] == 1.0 and p_bond['bond_type'] == 2.0:
                        if p_bond['stereo'] in (Chem.BondStereo.STEREOZ, Chem.BondStereo.STEREOE):
                            atom1, atom2 = bond_key
                            sym1 = reactant_atom_info.get(atom1, {}).get('symbol', 'X')
                            sym2 = reactant_atom_info.get(atom2, {}).get('symbol', 'X')

                            config = "Z (cis)" if p_bond['stereo'] == Chem.BondStereo.STEREOZ else "E (trans)"
                            change_desc = {
                                'bond': f"{sym1}:{atom1}={sym2}:{atom2}",
                                'type': 'E/Z',
                                'change': f"Single bond oxidized to {config} double bond"
                            }
                            explanation['stereochemistry_changes'].append(change_desc)
                            reaction_center.add(atom1)
                            reaction_center.add(atom2)

                    # Case 3: Stereo double bond -> Single/triple (loss of E/Z)
                    elif r_bond['bond_type'] == 2.0 and p_bond['bond_type'] in (1.0, 3.0):
                        if r_bond['stereo'] in (Chem.BondStereo.STEREOZ, Chem.BondStereo.STEREOE):
                            atom1, atom2 = bond_key
                            sym1 = reactant_atom_info.get(atom1, {}).get('symbol', 'X')
                            sym2 = reactant_atom_info.get(atom2, {}).get('symbol', 'X')

                            config = "Z (cis)" if r_bond['stereo'] == Chem.BondStereo.STEREOZ else "E (trans)"
                            new_bond = "single" if p_bond['bond_type'] == 1.0 else "triple"
                            change_desc = {
                                'bond': f"{sym1}:{atom1}-{sym2}:{atom2}",
                                'type': 'E/Z',
                                'change': f"{config} double bond changed to {new_bond} bond (loss of E/Z)"
                            }
                            explanation['stereochemistry_changes'].append(change_desc)
                            reaction_center.add(atom1)
                            reaction_center.add(atom2)

                    # Case 4: E/Z isomerization
                    elif r_bond['bond_type'] == 2.0 and p_bond['bond_type'] == 2.0:
                        if r_bond['stereo'] != p_bond['stereo']:
                            if (r_bond['stereo'] in (Chem.BondStereo.STEREOZ, Chem.BondStereo.STEREOE) and
                                    p_bond['stereo'] in (Chem.BondStereo.STEREOZ, Chem.BondStereo.STEREOE)):
                                atom1, atom2 = bond_key
                                sym1 = reactant_atom_info.get(atom1, {}).get('symbol', 'X')
                                sym2 = reactant_atom_info.get(atom2, {}).get('symbol', 'X')

                                r_config = "Z (cis)" if r_bond['stereo'] == Chem.BondStereo.STEREOZ else "E (trans)"
                                p_config = "Z (cis)" if p_bond['stereo'] == Chem.BondStereo.STEREOZ else "E (trans)"
                                change_desc = {
                                    'bond': f"{sym1}:{atom1}={sym2}:{atom2}",
                                    'type': 'E/Z',
                                    'change': f"Isomerization: {r_config} -> {p_config}"
                                }
                                explanation['stereochemistry_changes'].append(change_desc)
                                reaction_center.add(atom1)
                                reaction_center.add(atom2)

                # Case 5: New stereo double bond formed (not from existing bond)
                elif not r_bond and p_bond:
                    if p_bond['bond_type'] == 2.0 and p_bond['stereo'] in (
                    Chem.BondStereo.STEREOZ, Chem.BondStereo.STEREOE):
                        atom1, atom2 = bond_key
                        sym1 = product_atom_info.get(atom1, {}).get('symbol', 'X')
                        sym2 = product_atom_info.get(atom2, {}).get('symbol', 'X')

                        config = "Z (cis)" if p_bond['stereo'] == Chem.BondStereo.STEREOZ else "E (trans)"
                        change_desc = {
                            'bond': f"{sym1}:{atom1}={sym2}:{atom2}",
                            'type': 'E/Z',
                            'change': f"New {config} double bond formed"
                        }
                        explanation['stereochemistry_changes'].append(change_desc)
                        reaction_center.add(atom1)
                        reaction_center.add(atom2)

        explanation['reaction_center'] = reaction_center

        # Analyze bond changes using the transformation matrix
        if hasattr(classifier, 'transformation_matrix'):
            T = classifier.transformation_matrix
            lost_atoms = []
            added_atoms_list = []

            for i in range(len(classifier.be_matrix_reactants)):
                in_reactants = (np.any(classifier.be_matrix_reactants[i, :] != 0) or
                                np.any(classifier.be_matrix_reactants[:, i] != 0))

                in_products = (np.any(classifier.be_matrix_products[i, :] != 0) or
                               np.any(classifier.be_matrix_products[:, i] != 0))

                if in_reactants and not in_products:
                    lost_atoms.append(str(i))

                elif not in_reactants and in_products:
                    added_atoms_list.append(str(i))

            for i in range(len(T)):
                for j in range(i + 1, len(T[0])):
                    if T[i][j] != 0:
                        # Bond changed
                        atom_i = classifier.transformation_mapping[i]
                        atom_j = classifier.transformation_mapping[j]
                        r_idx_i = np.where(classifier.mapping_diagonal == atom_i)[0][0]
                        r_idx_j = np.where(classifier.mapping_diagonal == atom_j)[0][0]
                        sym_i = _element_symbol(classifier.atoms_diagonal[r_idx_i])
                        sym_j = _element_symbol(classifier.atoms_diagonal[r_idx_j])

                        if T[i][j] > 0:
                            reactant_connection = classifier.be_matrix_reactants[r_idx_i][r_idx_j]
                            product_connection = classifier.be_matrix_products[r_idx_i][r_idx_j]
                            if reactant_connection == 0:
                                bond_type = bond_name[product_connection]
                                explanation['bond_changes'].append(
                                    f"{sym_i}-{sym_j} {bond_type} bond formed between atoms {atom_i} and {atom_j}.")
                            else:
                                reactant_bond_type = bond_name[reactant_connection]
                                product_bond_type = bond_name[product_connection]
                                explanation['bond_changes'].append(
                                    f"{sym_i}-{sym_j} bond between atoms {atom_i} and {atom_j} oxidized "
                                    f"from {reactant_bond_type} to {product_bond_type}.")
                        else:
                            reactant_connection = classifier.be_matrix_reactants[r_idx_i][r_idx_j]
                            product_connection = classifier.be_matrix_products[r_idx_i][r_idx_j]
                            bond_type = bond_name[reactant_connection]
                            if product_connection == 0:
                                in_product_i = len(np.where(classifier.be_matrix_products[r_idx_i] != 0)[0])
                                in_product_j = len(np.where(classifier.be_matrix_products[r_idx_j] != 0)[0])

                                if in_product_i == 0 and in_product_j == 0:
                                    lost_atoms.append(str(r_idx_i))
                                    lost_atoms.append(str(r_idx_j))
                                    i_nbs = np.where(classifier.be_matrix_reactants[r_idx_i] != 0)
                                    i_no_nbs = True
                                    for nb in i_nbs:
                                        if len(np.where(classifier.be_matrix_products[nb] != 0)[0]) != 0:
                                            i_no_nbs = False
                                    j_nbs = np.where(classifier.be_matrix_reactants[r_idx_j] != 0)
                                    j_no_nbs = True
                                    for nb in j_nbs:
                                        if len(np.where(classifier.be_matrix_products[nb] != 0)[0]) != 0:
                                            j_no_nbs = False
                                    if i_no_nbs or j_no_nbs:
                                        explanation['leaving_groups'].append(
                                            f"{sym_i}-{sym_j} {bond_type} bond part of leaving group.")
                                    else:
                                        explanation['bond_changes'].append(
                                            f"{sym_i}-{sym_j} {bond_type} bond broken "
                                            f"between atoms {atom_i} and {atom_j}.")

                                else:
                                    if in_product_i == 0:
                                        lost_atoms.append(str(r_idx_i))
                                    if in_product_j == 0:
                                        lost_atoms.append(str(r_idx_j))
                                    explanation['bond_changes'].append(
                                        f"{sym_i}-{sym_j} {bond_type} bond broken between atoms {atom_i} and {atom_j}.")
                            else:
                                reactant_bond_type = bond_name[reactant_connection]
                                product_bond_type = bond_name[product_connection]
                                explanation['bond_changes'].append(
                                    f"{sym_i}-{sym_j} bond between atoms {atom_i} and {atom_j} reduced "
                                    f"from {reactant_bond_type} to {product_bond_type}.")
            if len(lost_atoms) > 0:
                lost_atoms = sorted(list(set(lost_atoms)))
                lost_symbols = [str(_element_symbol(classifier.atoms_diagonal[int(aidx)])) for aidx in lost_atoms]
                explanation['leaving_groups'].append(
                    f"Atom types {str(','.join(lost_symbols))}, with atom indices {','.join(lost_atoms)}")
            if len(added_atoms_list) > 0:
                added_atoms_list = sorted(list(set(added_atoms_list)))
                added_symbols = [str(_element_symbol(classifier.atoms_diagonal[int(aidx)])) for aidx in added_atoms_list]
                explanation['added_atoms'].append(
                    f"Atom types {str(','.join(added_symbols))}, with atom indices {','.join(added_atoms_list)}")

        # Get reaction classification info
        if hasattr(self, 'reaction_class'):
            explanation['classification']['class'] = self.reaction_class
        if hasattr(self, 'name'):
            explanation['classification']['name'] = self.name

        if len(self.reaction_info["FG_REACTANTS"]) > 0:
            explanation['functional_groups']['reacting'] = ", ".join(self.reaction_info["FG_REACTANTS"])
        if len(self.reaction_info["FG_PRODUCTS"]) > 0:
            explanation['functional_groups']['formed'] = ", ".join(self.reaction_info["FG_PRODUCTS"])

        from importlib.resources import files

        path = files(f"{__package__}.data").joinpath("named_rings.json")
        with path.open('r') as f:
            ring_dict = json.load(f)

        reactant_rings = set(self.reaction_info['PARTICIPATING_RINGS_REACTANTS'])
        product_rings = set(self.reaction_info['PARTICIPATING_RINGS_PRODUCTS'])

        # Rings that were formed (in products but not in reactants)
        rings_formed = list(product_rings - reactant_rings)

        # Rings that disappeared (in reactants but not in products)
        rings_broken = list(reactant_rings - product_rings)

        # Rings that were preserved (in both)
        rings_preserved = list(reactant_rings & product_rings)

        if len(rings_formed) > 0:
            explanation['rings']['formed'] = [ring_dict.get(r, r) for r in rings_formed]
        if len(rings_broken) > 0:
            explanation['rings']['broken'] = [ring_dict.get(r, r) for r in rings_broken]
        if len(rings_preserved) > 0:
            explanation['rings']['preserved'] = [ring_dict.get(r, r) for r in rings_preserved]

        # Generate template (RDChiral only)
        try:
            explanation['template'] = get_reaction_template(
                self.classifier.sanitized_mapped_reaction,
                radius_reactants=1,
                radius_products=0,
            )
        except Exception as e:
            explanation['template'] = f"Error generating template: {str(e)}"

        return explanation

