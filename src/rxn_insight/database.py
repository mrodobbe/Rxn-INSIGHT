import pandas as pd
from tqdm import tqdm
from rxn_insight.reaction import Reaction
from rxn_insight.utils import curate_smirks, get_reaction_template, make_rdkit_fp
from rxnmapper import RXNMapper
from typing import Union, Tuple, List
from importlib import resources


class Database:
    def __init__(self, df: Union[pd.DataFrame, None] = None):
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
        all_cols = [solvent_column, reagent_column, catalyst_column, yield_column, ref_column]
        df_protocol = df.rename(columns={reaction_column: "REACTION",
                                             solvent_column: "SOLVENT",
                                             reagent_column: "REAGENT",
                                             catalyst_column: "CATALYST",
                                             yield_column: "YIELD",
                                             ref_column: "REF"})
        for col in all_cols:
            if col not in df_protocol.keys():
                df_protocol[col] = ""
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
        df_csv = pd.read_csv(fname, index_col=None)
        all_cols = [solvent_column, reagent_column, catalyst_column, yield_column, ref_column]
        df_protocol = df_csv.rename(columns={reaction_column: "REACTION",
                                             solvent_column: "SOLVENT",
                                             reagent_column: "REAGENT",
                                             catalyst_column: "CATALYST",
                                             yield_column: "YIELD",
                                             ref_column: "REF"})
        for col in all_cols:
            if col not in df_protocol.keys():
                df_protocol[col] = ""
        df, skipped_reactions = analyze_reactions(df_protocol)
        self.df = df
        self.skipped_reactions = skipped_reactions
        self.class_distribution = calculate_class_distribution(df)
        self.name_distribution = calculate_name_distribution(df)

        return df

    def save_to_parquet(self, fname: str):
        self.df.to_parquet(f"{fname}.gzip")

    def save_to_csv(self, fname: str):
        self.df.to_csv(f"{fname}.gzip")

    def save_to_excel(self, fname: str):
        self.df.to_excel(f"{fname}.xlsx")

    def get_class_distribution(self):
        return self.class_distribution

    def get_name_distribution(self):
        return self.name_distribution


def analyze_reactions(df: pd.DataFrame) -> Tuple[pd.DataFrame, List[str]]:
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
    names_dict = {"NAME": [], "COUNT": []}
    all_names = df["NAME"].to_list()
    unique_names = list(set(all_names))
    for name in unique_names:
        names_dict["NAME"].append(name)
        names_dict["COUNT"].append(all_names.count(name))

    df_names = pd.DataFrame(names_dict)

    return df_names
