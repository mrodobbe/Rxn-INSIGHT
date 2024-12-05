import pandas as pd
from tqdm import tqdm
from rxn_insight.reaction import Reaction
from rxn_insight.utils import curate_smirks, get_reaction_template, make_rdkit_fp
from rxnmapper import RXNMapper
from typing import Union, Tuple, List
from importlib import resources


class Database:

    """
    A class to manage and analyze reaction datasets, providing functionalities
    for creating databases, analyzing reactions, and saving results.

    Example:
        >>> from rxn_insight.database import Database
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
        >>> db = Database()
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

        if df is None:
            self.df = pd.DataFrame({})
        else:
            self.df = df
        self.skipped_reactions = []
        self.class_distribution = pd.DataFrame({})
        self.name_distribution = pd.DataFrame({})

    def create_database_from_df(
            self,
            df: pd.DataFrame,
            reaction_column: str,
            solvent_column: str = "SOLVENT",
            reagent_column: str = "REAGENT",
            catalyst_column: str = "CATALYST",
            yield_column: str = "YIELD",
            ref_column: str = "REF"
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

        df, skipped_reactions = analyze_reactions(df_protocol)
        self.df = df
        self.skipped_reactions = skipped_reactions
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
            ref_column: str = "REF"
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

        df, skipped_reactions = analyze_reactions(df_protocol)
        self.df = df
        self.skipped_reactions = skipped_reactions
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

        return self.class_distribution

    def get_name_distribution(self):

        """
        Retrieves the distribution of reaction names in the database.

        Returns:
            A DataFrame summarizing the reaction name distribution.
        """

        return self.name_distribution


def analyze_reactions(df: pd.DataFrame) -> Tuple[pd.DataFrame, List[str]]:

    """
    Analyzes a DataFrame of reactions to extract detailed information.

    Args:
        df: A DataFrame with reaction data.

    Returns:
        A tuple containing the updated DataFrame and a list of skipped reactions.

    """

    headers = [
        'REACTANTS', 'PRODUCTS', 'SANITIZED_REACTION', 'MAPPED_REACTION', 'N_REACTANTS', 'N_PRODUCTS',
        'FG_REACTANTS', 'FG_PRODUCTS', 'PARTICIPATING_RINGS_REACTANTS', 'PARTICIPATING_RINGS_PRODUCTS',
        'ALL_RINGS_PRODUCTS', 'BY-PRODUCTS', 'CLASS', 'NAME',
        'TAG', 'TAG2', 'SCAFFOLD', 'rxn_str_patt_fp', 'rxn_dif_patt_fp',
        'rxn_str_morgan_fp', 'rxn_dif_morgan_fp', 'TEMPLATE'
    ]
    df[headers] = ""

    rxn_mapper = RXNMapper()
    with resources.path(
            f"{__package__}.json", "functional_groups.json"
    ) as path:
        fg_db = pd.read_json(path, orient="records", lines=True)

    with resources.path(
            f"{__package__}.json", "smirks.json"
    ) as path:
        smirks = pd.read_json(path, orient="records", lines=True)

    smirks_db = curate_smirks(smirks)
    bad_reactions = []

    for i in tqdm(df.index, desc="Analyzing all reactions..."):
        try:
            rxn = Reaction(df["REACTION"][i],
                           solvent=df["SOLVENT"][i],
                           reagent=df["REAGENT"][i],
                           catalyst=df["CATALYST"][i], ref=df["REF"][i], rxn_mapper=rxn_mapper,
                           smirks=smirks_db, fg=fg_db)
            ri = rxn.get_reaction_info()
            for header in headers:
                if header in list(ri.keys()):
                    if type(ri[header]) is list:
                        df.loc[i, header] = ".".join(ri[header])
                    else:
                        df.loc[i, header] = ri[header]

            df.loc[i, "SANITIZED_REACTION"] = ri["REACTION"]
            df.loc[i, "TAG2"] = rxn.give_broad_tag()
            df.loc[i, "TEMPLATE"] = get_reaction_template(ri["MAPPED_REACTION"], 2, 2)
            df.loc[i, "REACTANTS"], df.loc[i, "PRODUCTS"] = ri["REACTION"].split(">>")

            df.loc[i, "rxn_str_patt_fp"] = make_rdkit_fp(ri["REACTION"], fp="MACCS", concatenate=True)
            df.loc[i, "rxn_dif_patt_fp"] = make_rdkit_fp(ri["REACTION"], fp="MACCS", concatenate=False)
            df.loc[i, "rxn_str_morgan_fp"] = make_rdkit_fp(ri["REACTION"], fp="Morgan", concatenate=True)
            df.loc[i, "rxn_dif_morgan_fp"] = make_rdkit_fp(ri["REACTION"], fp="Morgan", concatenate=False)

        except KeyboardInterrupt:
            raise
        except Exception as e:
            print(f"Error at reaction nr {i}, {df['REACTION'][i]}", e)
            bad_reactions.append(df['REACTION'][i])
            df = df.drop(i)

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
