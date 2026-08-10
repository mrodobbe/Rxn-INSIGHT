from __future__ import annotations

import numpy as np
import pandas as pd
from tqdm import tqdm
from rxn_insight.reaction import Reaction
from rxn_insight.classification import ReactionClassifier
from rxn_insight.utils import curate_smirks, get_reaction_template, make_rdkit_fp
from typing import Union, Tuple, List, Dict, Optional
from importlib import resources
import multiprocessing as mp
from joblib import Parallel, delayed


class Database:
    """
    A class to manage and analyze reaction datasets, providing functionalities
    for creating databases, analyzing reactions, and saving results.

    Example:
        >>> import rxn_insight as ri
        >>> import pandas as pd
        >>> # Create a sample DataFrame
        >>> data = {
        ...     "reaction": ["OB(O)c1ccccc1.Brc1ccccc1>>c1ccc(-c2ccccc2)cc1"],
        ...     "solvent": ["CN(C)C=O"],
        ...     "reagent": ["F[Cs]"],
        ...     "catalyst": ["[Pd]"],
        ...     "yield": [85],
        ...     "reference": ["Ref1"]
        ... }
        >>> df = pd.DataFrame(data)
        >>> # Initialize a Database object
        >>> db = ri.Database()
        >>> # Create a database from the DataFrame
        >>> reaction_df = db.create_database_from_df(
        ...     df,
        ...     reaction_column="reaction",
        ...     solvent_column="solvent",
        ...     reagent_column="reagent",
        ...     catalyst_column="catalyst",
        ...     yield_column="yield",
        ...     ref_column="reference"
        ...     )

    """

    def __init__(self, df: Union[pd.DataFrame, None] = None):

        """
        Initializes a Database object with an optional DataFrame.

        Args:
            df: An optional pandas DataFrame containing reaction data.

        """

        self.df = pd.DataFrame({})
        self.class_distribution = pd.DataFrame({})
        self.name_distribution = pd.DataFrame({})

        if df is not None:
            self.df = df
            self.class_distribution = self.get_class_distribution()
            self.name_distribution = self.get_name_distribution()

        self.skipped_reactions = []

    def create_database_from_df(
            self,
            df: pd.DataFrame,
            reaction_column: str,
            solvent_column: str = "SOLVENT",
            reagent_column: str = "REAGENT",
            catalyst_column: str = "CATALYST",
            yield_column: str = "YIELD",
            ref_column: str = "REF",
            classify: bool = True,
            add_fp: bool = True,
            n_jobs: Optional[int] = None,
    ) -> pd.DataFrame:

        """
        Creates a reaction database from a given DataFrame.

        Args:
            df: A DataFrame containing reaction data.
            reaction_column: Name of the column containing reaction SMILES.
            solvent_column: Name of the solvent column (default: "SOLVENT").
            reagent_column: Name of the reagent column (default: "REAGENT").
            catalyst_column: Name of the catalyst column (default: "CATALYST").
            yield_column: Name of the yield column (default: "YIELD").
            ref_column: Name of the reference column (default: "REF").

        Returns:
            A DataFrame with analyzed reaction data.
        """

        all_cols = ["SOLVENT", "REAGENT", "CATALYST", "YIELD", "REF"]
        df_protocol = df.rename(columns={reaction_column: "REACTION",
                                         solvent_column: "SOLVENT",
                                         reagent_column: "REAGENT",
                                         catalyst_column: "CATALYST",
                                         yield_column: "YIELD",
                                         ref_column: "REF"})
        for col in all_cols:
            if col not in df_protocol.keys():
                df_protocol[col] = "not-reported"

        df, skipped_reactions = analyze_reactions(df_protocol, classify=classify, add_fp=add_fp, n_jobs=n_jobs)
        self.df = df
        self.skipped_reactions = skipped_reactions
        if classify:
            self.class_distribution = calculate_class_distribution(df)
        self.name_distribution = calculate_name_distribution(df)

        return df

    def create_database_from_csv(
            self,
            fname: str,
            reaction_column: str,
            solvent_column: str = "SOLVENT",
            reagent_column: str = "REAGENT",
            catalyst_column: str = "CATALYST",
            yield_column: str = "YIELD",
            ref_column: str = "REF",
            classify: bool = True,
            add_fp: bool = True,
            n_jobs: Optional[int] = None,
    ) -> pd.DataFrame:

        """
        Creates a reaction database from a CSV file.

        Args:
            fname: Path to the CSV file.
            reaction_column: Name of the column containing reaction SMILES.
            solvent_column: Name of the solvent column (default: "SOLVENT").
            reagent_column: Name of the reagent column (default: "REAGENT").
            catalyst_column: Name of the catalyst column (default: "CATALYST").
            yield_column: Name of the yield column (default: "YIELD").
            ref_column: Name of the reference column (default: "REF").

        Returns:
            A DataFrame with analyzed reaction data.
        """

        df_csv = pd.read_csv(fname, index_col=None)
        all_cols = ["SOLVENT", "REAGENT", "CATALYST", "YIELD", "REF"]
        df_protocol = df_csv.rename(columns={reaction_column: "REACTION",
                                             solvent_column: "SOLVENT",
                                             reagent_column: "REAGENT",
                                             catalyst_column: "CATALYST",
                                             yield_column: "YIELD",
                                             ref_column: "REF"})
        for col in all_cols:
            if col not in df_protocol.keys():
                df_protocol[col] = "not-reported"

        df, skipped_reactions = analyze_reactions(df_protocol, classify=classify, add_fp=add_fp, n_jobs=n_jobs)
        self.df = df
        self.skipped_reactions = skipped_reactions
        if classify:
            self.class_distribution = calculate_class_distribution(df)
        self.name_distribution = calculate_name_distribution(df)

        return df

    def save_to_parquet(self, fname: str):

        """
        Saves the reaction database to a Parquet file.

        Args:
            fname: The name of the output file (without extension).
        """

        self.df.to_parquet(f"{fname}.gzip")

    def save_to_csv(self, fname: str):

        """
        Saves the reaction database to a CSV file.

        Args:
            fname: The name of the output file (without extension).
        """

        self.df.to_csv(f"{fname}.csv")

    def save_to_excel(self, fname: str):

        """
        Saves the reaction database to an Excel file.

        Args:
            fname: The name of the output file (without extension).
        """

        self.df.to_excel(f"{fname}.xlsx")

    def get_class_distribution(self):

        """
        Retrieves the class distribution of reactions in the database.

        Returns:
            A DataFrame summarizing the reaction class distribution.
        """
        if len(self.class_distribution) == 0:
            self.class_distribution = calculate_class_distribution(self.df)
        return self.class_distribution

    def get_name_distribution(self):

        """
        Retrieves the distribution of reaction names in the database.

        Returns:
            A DataFrame summarizing the reaction name distribution.
        """
        if len(self.name_distribution) == 0:
            self.name_distribution = calculate_name_distribution(self.df)

        return self.name_distribution


def process_single_reaction(
        idx: int,
        reaction: str,
        solvent: str,
        reagent: str,
        catalyst: str,
        ref: str,
        rxn_mapper,
        smirks_db,
        fg_db,
        classify: bool = True,
        add_fp: bool = True,
) -> Tuple[int, Optional[Dict], Optional[str]]:
    """
    Process a single reaction and return results or error.

    Returns:
        Tuple of (index, result_dict, error_reaction)
    """
    try:
        rxn = Reaction(
            reaction,
            solvent=solvent,
            reagent=reagent,
            catalyst=catalyst,
            ref=ref,
            rxn_mapper=rxn_mapper,
            smirks=smirks_db,
            fg=fg_db,
            classify=classify,
        )
        if classify:
            rinfo = rxn.get_reaction_info()

            # Prepare result dictionary
            result = {}
            headers = [
                'REACTANTS', 'PRODUCTS', 'SANITIZED_REACTION', 'MAPPED_REACTION', 'N_REACTANTS', 'N_PRODUCTS',
                'FG_REACTANTS', 'FG_PRODUCTS', 'PARTICIPATING_RINGS_REACTANTS', 'PARTICIPATING_RINGS_PRODUCTS',
                'ALL_RINGS_PRODUCTS', 'BY-PRODUCTS', 'CLASS', 'NAME',
                'TAG', 'TAG2', 'SCAFFOLD', 'rxn_str_patt_fp', 'rxn_dif_patt_fp',
                'rxn_str_morgan_fp', 'rxn_dif_morgan_fp', 'TEMPLATE',
            ]

            for header in headers:
                if header in rinfo:
                    if isinstance(rinfo[header], list):
                        result[header] = ".".join(rinfo[header])
                    else:
                        result[header] = rinfo[header]

            # Additional processing
            result["SANITIZED_REACTION"] = rinfo["REACTION"]
            result["TAG2"] = rxn.give_broad_tag()
            result["TEMPLATE"] = get_reaction_template(rinfo["MAPPED_REACTION"], 2, 2)

            result["REACTANTS"], result["PRODUCTS"] = rinfo["REACTION"].split(">>")

            if add_fp:
                result["rxn_str_patt_fp"] = make_rdkit_fp(rinfo["REACTION"], fp="MACCS", concatenate=True)
                result["rxn_dif_patt_fp"] = make_rdkit_fp(rinfo["REACTION"], fp="MACCS", concatenate=False)
                result["rxn_str_morgan_fp"] = make_rdkit_fp(rinfo["REACTION"], fp="Morgan", concatenate=True)
                result["rxn_dif_morgan_fp"] = make_rdkit_fp(rinfo["REACTION"], fp="Morgan", concatenate=False)
        else:
            result = {"NAME": rxn.get_name()}

            if add_fp:
                result["rxn_str_patt_fp"] = make_rdkit_fp(rxn.reaction, fp="MACCS", concatenate=True)
                result["rxn_dif_patt_fp"] = make_rdkit_fp(rxn.reaction, fp="MACCS", concatenate=False)
                result["rxn_str_morgan_fp"] = make_rdkit_fp(rxn.reaction, fp="Morgan", concatenate=True)
                result["rxn_dif_morgan_fp"] = make_rdkit_fp(rxn.reaction, fp="Morgan", concatenate=False)

        return idx, result, None

    except Exception as e:
        print(f"Error at reaction nr {idx}, {reaction}", e)
        return idx, None, reaction


def analyze_reactions(
        df: pd.DataFrame,
        classify: bool = True,
        add_fp: bool = True,
        n_jobs: Optional[int] = None,
) -> Tuple[pd.DataFrame, List[str]]:
    """
    Analyzes a DataFrame of reactions to extract detailed information.

    Args:
        df: A DataFrame with reaction data.

    Returns:
        A tuple containing the updated DataFrame and a list of skipped reactions.

    """

    if n_jobs is None:
        n_jobs = max(1, int(0.5 * mp.cpu_count()))

    headers = [
        'REACTANTS', 'PRODUCTS', 'SANITIZED_REACTION', 'MAPPED_REACTION', 'N_REACTANTS', 'N_PRODUCTS',
        'FG_REACTANTS', 'FG_PRODUCTS', 'PARTICIPATING_RINGS_REACTANTS', 'PARTICIPATING_RINGS_PRODUCTS',
        'ALL_RINGS_PRODUCTS', 'BY-PRODUCTS', 'CLASS', 'NAME',
        'TAG', 'TAG2', 'SCAFFOLD', 'rxn_str_patt_fp', 'rxn_dif_patt_fp',
        'rxn_str_morgan_fp', 'rxn_dif_morgan_fp', 'TEMPLATE',
    ]
    for col in headers:
        df[col] = pd.array([None] * len(df), dtype=object)

    from rxnmapper import RXNMapper
    rxn_mapper = RXNMapper()
    fg_db = pd.read_json(
        resources.files(f"{__package__}.data").joinpath("functional_groups.json"),
        orient="records", lines=True,
    )
    smirks = pd.read_json(
        resources.files(f"{__package__}.data").joinpath("smirks.json"),
        orient="records", lines=True,
    )

    smirks_db = curate_smirks(smirks)
    bad_reactions = []

    reaction_data = [
        (
            idx,
            df.loc[idx, "REACTION"],
            df.loc[idx, "SOLVENT"],
            df.loc[idx, "REAGENT"],
            df.loc[idx, "CATALYST"],
            df.loc[idx, "REF"]
        )
        for idx in df.index
    ]

    results = Parallel(n_jobs=n_jobs)(
        delayed(process_single_reaction)(
            idx, reaction, solvent, reagent, catalyst, ref,
            rxn_mapper, smirks_db, fg_db, classify, add_fp,
        )
        for idx, reaction, solvent, reagent, catalyst, ref in tqdm(
            reaction_data, desc=f"Analyzing reactions (using {n_jobs} cores)..."
        )
    )

    indices_to_drop = []
    for idx, result, error_reaction in results:
        if result is not None:
            for key, value in result.items():
                df.loc[idx, key] = value
        else:
            indices_to_drop.append(idx)
            bad_reactions.append(error_reaction)

    df = df.drop(indices_to_drop)

    return df, bad_reactions


def calculate_class_distribution(df: pd.DataFrame) -> pd.DataFrame:
    """
    Calculates the distribution of reaction classes.

    Args:
        df: A DataFrame containing reaction data.

    Returns:
        A DataFrame summarizing reaction class counts.
    """

    class_dict = {"CLASS": [], "COUNT": []}
    classes = ['Acylation',
               'Heteroatom Alkylation and Arylation',
               'Aromatic Heterocycle Formation',
               'C-C Coupling',
               'Deprotection',
               'Protection',
               'Functional Group Interconversion',
               'Functional Group Addition',
               'Reduction',
               'Oxidation',
               'Miscellaneous']
    all_classes = df["CLASS"].to_list()
    for c in classes:
        class_dict["CLASS"].append(c)
        class_dict["COUNT"].append(all_classes.count(c))

    df_class = pd.DataFrame(class_dict)

    return df_class


def _rdchiral_templates_for_chunk(
        chunk: List[str],
        radii: List[Tuple[int, int]],
) -> List[Dict]:
    """Worker: extract RDChiral templates for a chunk of reactions."""
    results = []
    for mapped_reaction in chunk:
        row: Dict = {}
        for (r_r, r_p) in radii:
            try:
                row[(r_r, r_p)] = get_reaction_template(
                    mapped_reaction, radius_reactants=r_r, radius_products=r_p
                )
            except Exception:
                row[(r_r, r_p)] = None
        results.append(row)
    return results


def extract_rdchiral_templates_batch(
        mapped_reactions: Union[List[str], "pd.Series"],
        radii: List[Union[int, Tuple[int, int]]] = None,
        n_jobs: Optional[int] = None,
        chunk_size: Optional[int] = None,
        progress: bool = True,
) -> pd.DataFrame:
    """Extract RDChiral reaction templates for a large number of pre-mapped reactions.

    A parallelised batch wrapper around :func:`get_reaction_template` (which
    calls ``rdchiral.template_extractor.extract_from_reaction`` internally).
    Only the radius parameters vary; unlike the Rxn-INSIGHT extractor there is
    no ``ring_info`` toggle because RDChiral always encodes ring membership in
    its SMARTS atom specs.

    Args:
        mapped_reactions: Iterable of **atom-mapped** reaction SMILES strings.
            RDChiral requires atom maps — unmapped reactions will yield ``None``.
        radii: Radii to compute templates for (default: ``[0, 1, 2, 3]``).
            Each element may be an ``int`` (symmetric: same radius on both
            sides) or a ``(r_reactants, r_products)`` tuple for asymmetric
            templates.
        n_jobs: Parallel workers. ``None`` → half the available CPU cores;
            ``-1`` → all cores.
        chunk_size: Reactions per worker call. ``None`` → chosen automatically
            as ``max(50, total // (n_jobs * 8))``.
        progress: Show a tqdm progress bar (default: True).

    Returns:
        DataFrame with columns:

        * ``MAPPED_REACTION`` — the input SMILES.
        * ``RDCHIRAL_rr{r_r}rp{r_p}`` — template string for each radius pair,
          or ``None`` if extraction failed.

    Example:
        >>> reactions = [
        ...     "Cl[C:2]([CH3:1])=[O:3].[NH2:4][c:5]1ccccc1"
        ...     ">>[CH3:1][C:2](=[O:3])[NH:4][c:5]1ccccc1",
        ... ]
        >>> df = extract_rdchiral_templates_batch(reactions, radii=[0, 1, 2])
        >>> df.columns.tolist()
        ['MAPPED_REACTION', 'RDCHIRAL_rr0rp0', 'RDCHIRAL_rr1rp1', 'RDCHIRAL_rr2rp2']
    """
    if radii is None:
        radii = [0, 1, 2, 3]
    if n_jobs is None:
        n_jobs = max(1, mp.cpu_count() // 2)

    # Normalize to (r_r, r_p) tuples, deduplicate while preserving order
    normalized: List[Tuple[int, int]] = [(r, r) if isinstance(r, int) else (r[0], r[1]) for r in radii]
    seen_pairs: set = set()
    normalized_ordered: List[Tuple[int, int]] = []
    for pair in normalized:
        if pair not in seen_pairs:
            seen_pairs.add(pair)
            normalized_ordered.append(pair)

    reactions_list = list(mapped_reactions)
    n = len(reactions_list)

    effective_jobs = n_jobs if n_jobs > 0 else mp.cpu_count()
    if chunk_size is None:
        chunk_size = max(50, n // (effective_jobs * 8))

    chunks = [reactions_list[i: i + chunk_size] for i in range(0, n, chunk_size)]

    raw_chunks: List[List[Dict]] = Parallel(n_jobs=n_jobs)(
        delayed(_rdchiral_templates_for_chunk)(chunk, normalized_ordered)
        for chunk in tqdm(chunks, desc="Extracting RDChiral templates", disable=not progress)
    )

    raw: List[Dict] = [item for chunk_res in raw_chunks for item in chunk_res]

    col_names = [f"RDCHIRAL_rr{r_r}rp{r_p}" for (r_r, r_p) in normalized_ordered]
    data: Dict[str, list] = {"MAPPED_REACTION": np.array(reactions_list, dtype=object)}
    for (r_r, r_p), col in zip(normalized_ordered, col_names):
        data[col] = np.array([r.get((r_r, r_p)) for r in raw], dtype=object)

    return pd.DataFrame(data)


def calculate_name_distribution(df: pd.DataFrame) -> pd.DataFrame:
    """
    Calculates the distribution of reaction names.

    Args:
        df: A DataFrame containing reaction data.

    Returns:
        A DataFrame summarizing reaction name counts.
    """

    names_dict = {"NAME": [], "COUNT": []}
    all_names = df["NAME"].to_list()
    unique_names = list(set(all_names))
    for name in unique_names:
        names_dict["NAME"].append(name)
        names_dict["COUNT"].append(all_names.count(name))

    df_names = pd.DataFrame(names_dict)

    return df_names
