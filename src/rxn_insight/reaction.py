"""Reaction module"""

import hashlib
import warnings
from typing import Optional

import numpy as np
import pandas as pd
from rdkit import Chem
from rxnmapper import RXNMapper
from tqdm import tqdm

from rxn_insight.classification import ReactionClassifier
from rxn_insight.utils import (
    atom_remover,
    curate_smirks,
    get_catalyst_ranking,
    get_fp,
    get_reagent_ranking,
    get_ring_systems,
    get_scaffold,
    get_similarity,
    get_solvent_ranking,
    maccs_fp,
    morgan_fp,
    remove_atom_mapping,
    sanitize_ring,
)


class Reaction:
    """This class handles operations related to chemical reactions."""

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
    ):
        """Initializes a Reaction object with reaction details and optional components.

        Args:
            reaction: A string representing the reaction in SMILES format.
            solvent: A string of solvents used in the reaction.
            reagent: A string of reagents used in the reaction.
            catalyst: A string of catalysts used in the reaction.
            ref: Reference or note associated with the reaction.
            rxn_mapper: An optional RXNMapper object for reaction mapping.
            keep_mapping: A boolean indicating if atom mappings should be preserved.
            smirks: A DataFrame containing SMIRKS transformations.
            fg: A DataFrame containing functional groups data.
        """
        self.reaction = ""
        self.solvent = solvent
        self.reagent = reagent
        self.catalyst = catalyst
        self.reference = ref
        self.read_reaction(reaction)
        if ":" in self.reaction and not keep_mapping:
            self.reaction = remove_atom_mapping(
                self.reaction
            )  # Remove atom mapping for consistency
        else:
            self.reaction = self.reaction
        self.smirks_db = smirks
        self.fg_db = fg
        self.classifier = ReactionClassifier(
            reaction, rxn_mapper=rxn_mapper, keep_mapping=keep_mapping
        )
        self.add_agents()
        self.reactants, self.products = self.classifier.sanitized_reaction.split(">>")
        self.mapped_reaction = self.classifier.sanitized_mapped_reaction
        self.reaction_class = ""
        self.template = self.classifier.template
        self.reaction_info: dict[str, tuple[str, ...] | str] = dict()
        self.tag = ""
        self.name = ""
        self.byproducts: tuple[str, ...] = tuple()
        self.scaffold = self.get_scaffold()
        self.neighbors = None
        self.suggested_solvent = ""
        self.suggested_catalyst = ""
        self.suggested_reagent = ""

    def read_reaction(self, reaction: str) -> None:
        """Processes the given reaction string to extract and format components.

        Args:
            reaction: The reaction string in SMILES format.
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
        """Adds agents from the classifier to the reagent list."""
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
            from importlib import resources

            with resources.path(
                f"{__package__}.json", "functional_groups.json"
            ) as path:
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
            from importlib import resources

            with resources.path(f"{__package__}.json", "smirks.json") as path:
                self.smirks_db = curate_smirks(
                    pd.read_json(path, orient="records", lines=True)
                )
        self.name = self.classifier.name_reaction(self.smirks_db)
        return self.name

    def get_reaction_info(self) -> dict[str, list[str] | str]:
        """Compiles detailed information about the reaction and its components."""
        if self.fg_db is None:
            from importlib import resources

            with resources.path(
                f"{__package__}.json", "functional_groups.json"
            ) as path:
                self.fg_db = pd.read_json(path, orient="records", lines=True)

        info_dict = self.classifier.get_reaction_center_info(self.fg_db)
        self.tag = info_dict["TAG"]
        self.reaction_class = info_dict["CLASS"]

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
            print("No similar reactions found...")
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
        df_tag["SOLVENT"].fillna("", inplace=True)
        df_tag["CATALYST"].fillna("", inplace=True)
        df_tag["REAGENT"].fillna("", inplace=True)
        max_similarity = df_tag["SIMILARITY"].max()
        df_tag = df_tag[df_tag["SIMILARITY"] > threshold].copy()
        print(
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


class Molecule:
    """This class reads in SMILES."""

    def __init__(self, smi: str):
        """Initializes a Molecule object with the SMILES string of a molecule.

        Args:
            smi: A string containing the SMILES representation of the molecule.
        """
        self.mol = Chem.MolFromSmiles(smi)
        self.smiles = Chem.MolToSmiles(self.mol)
        self.inchi = Chem.MolToInchi(self.mol)
        self.inchikey = Chem.MolToInchiKey(self.mol)
        self.functional_groups = None
        # self.rings = tuple() # Seems to be unused
        self.scaffold = get_scaffold(self.mol)
        self.maccs_fp = maccs_fp(self.mol)
        self.morgan_fp = morgan_fp(self.mol)
        self.reactions = None

    def search_reactions(self, df: pd.DataFrame) -> pd.DataFrame:
        """Searches for reactions involving the molecule as a product.

        Args:
            df: The DataFrame to search for reactions.
        """
        if "PRODUCT" in df.keys():
            dfc = df[df["PRODUCT"] == self.inchikey].copy()
        else:
            df["PRODUCT"] = ""
            for i in tqdm(df.index):
                try:
                    prod = df["REACTION"][i].split(">>")[1]
                    mol = Chem.MolFromSmiles(prod)
                    df["PRODUCT"][i] = Chem.MolToInchiKey(mol)
                except KeyboardInterrupt:
                    raise
                except Exception as e:
                    print(e)
                    continue
            dfc = df[df["PRODUCT"] == self.inchikey].copy()

        if "rxn_str_patt_fp" in dfc.keys():
            dfc = dfc.drop(columns=["rxn_str_patt_fp"])
        if "rxn_dif_patt_fp" in dfc.keys():
            dfc = dfc.drop(columns=["rxn_dif_patt_fp"])
        if "rxn_str_morgan_fp" in dfc.keys():
            dfc = dfc.drop(columns=["rxn_str_morgan_fp"])
        if "rxn_dif_morgan_fp" in dfc.keys():
            dfc = dfc.drop(columns=["rxn_dif_morgan_fp"])
        if "TAG" in dfc.keys():
            dfc = dfc.drop(columns=["TAG"])
        if "TAG2" in dfc.keys():
            dfc = dfc.drop(columns=["TAG2"])

        self.reactions = dfc

        return dfc

    def search_reactions_by_scaffold(
        self,
        df: pd.DataFrame,
        threshold: float = 0.5,
        max_return: int = 100,
        fp: str = "MACCS",
    ) -> pd.DataFrame:
        """Searches for reactions based on scaffold similarity.

        Args:
            df: DataFrame containing reactions to search.
            threshold: Similarity threshold to apply.
            max_return: Maximum number of reactions to return.
            fp: Type of fingerprint to use for similarity calculation.
        """
        dfc = df[df["SCAFFOLD"] == self.scaffold].copy()
        if len(dfc.index) == 0:
            print("No products with the same scaffold found!")
            return None

        if fp.lower() == "maccs":
            fps = [
                maccs_fp(Chem.MolFromSmiles(r.split(">>")[1]))
                for r in tqdm(dfc["REACTION"].tolist(), desc="Making fingerprints...")
            ]
            dfc["SIMILARITY"] = [get_similarity(self.maccs_fp, fp) for fp in fps]
        elif fp.lower() == "morgan":
            fps = [
                morgan_fp(Chem.MolFromSmiles(r.split(">>")[1]))
                for r in tqdm(dfc["REACTION"].tolist(), desc="Making fingerprints...")
            ]
            dfc["SIMILARITY"] = [get_similarity(self.morgan_fp, fp) for fp in fps]
        else:
            raise KeyError(
                f"Fingerprint choice {fp} is not supported. Select MACCS or Morgan."
            )

        df_tag = dfc.sort_values(by="SIMILARITY", ascending=False).copy()
        df_tag["SOLVENT"].fillna("", inplace=True)
        df_tag["CATALYST"].fillna("", inplace=True)
        df_tag["REAGENT"].fillna("", inplace=True)
        max_similarity = df_tag["SIMILARITY"].max()
        df_tag = df_tag[df_tag["SIMILARITY"] > threshold].copy()
        print(
            f"Product found with similarity of {max_similarity:.3f}. This will be our best match."
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

        return df_return

    def get_functional_groups(self, df: pd.DataFrame = None) -> list[str]:
        """Identifies and returns the functional groups present in the molecule.

        Args:
            df: Optional DataFrame containing functional group patterns; loads default if not provided.
        """
        if df is None:
            from importlib import resources

            with resources.path(
                f"{__package__}.json", "functional_groups.json"
            ) as path:
                df = pd.read_json(path, orient="records", lines=True)

        mol = self.mol
        atom_indices = np.array([atom.GetIdx() for atom in mol.GetAtoms()])
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

    def get_rings(self) -> list[str]:
        """Identifies and returns rings in the molecule."""
        mol = self.mol
        try:
            rs = get_ring_systems(mol, include_spiro=True)
        except:
            return []

        found_rings = []
        if len(rs) > 0:
            for k in range(len(rs)):
                found_rings.append(sanitize_ring(atom_remover(mol, [rs[k]])))
            return found_rings
        else:
            return []
