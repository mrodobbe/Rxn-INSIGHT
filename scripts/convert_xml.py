"""Script to convert USPTO data."""

import os
import xml.etree.ElementTree as ET

import pandas as pd


class Reaction:
    """Reaction class."""

    def __init__(self, tree):
        """Initialize the Reaction class with various properties extracted from a tree.

        Args:
            tree (xml.etree.ElementTree.Element): The XML element containing the
            reaction data.
        """
        self.tree = tree
        self.reference = self.get_reference()
        self.complete_reaction = self.get_rxn_smiles()
        self.reactants, self.products, self.reaction = self.get_reaction()
        self.reagents = self.get_reagents()
        self.solvents = self.reagents["solvent"]
        self.catalysts = self.reagents["catalyst"]
        self.paragraph = self.get_paragraphs()
        self.title = self.get_titles()

    def get_reference(self):
        """Extracts the reference from the XML tree.

        Returns:
        -------
            str: The reference text extracted from the XML.
        """
        return self.tree[0].getchildren()[0].text

    def get_rxn_smiles(self):
        """Extracts the SMILES representation of the reaction.

        Returns:
        -------
            str: The SMILES string representing the complete reaction.
        """
        return self.tree[1].text

    def get_titles(self):
        """Extracts the title from the XML tree.

        Returns:
        -------
            str: The title text extracted from the XML.
        """
        return self.tree[0][1].text

    def get_paragraphs(self):
        """Extracts paragraphs from the XML tree.

        Returns:
        -------
            str: The text of the last paragraph from the XML tree.
        """
        return self.tree[0][-1].text

    def get_reaction(self):
        """Processes the reaction SMILES string to extract reactants, products, the reaction.

        Returns:
        -------
            tuple: A tuple containing the reactants, products, and the formatted reaction string.
        """
        comps = self.complete_reaction.split(">")
        return comps[0], comps[-1], f"{comps[0]}>>{comps[-1]}"

    def get_reagents(self):
        """Extracts reagents, solvents, and catalysts from the XML tree.

        Returns:
        -------
            dict: A dictionary containing solvents and catalysts as keys and their respective SMILES as values.
        """
        reagents = self.tree[4].getchildren()
        if len(reagents) == 0:
            return {"solvent": "", "catalyst": ""}
        else:
            reagent_dict = dict()
            for reagent in reagents:
                reagent_type = reagent.attrib["role"]
                if reagent_type in reagent_dict:
                    for val in reagent.getchildren():
                        if "identifier" in val.tag:
                            reag_smile = val.attrib["value"]
                            break
                    try:
                        reagent_dict[reagent_type] = (
                            f"{reagent_dict[reagent_type]}.{reag_smile}"
                        )
                    except UnboundLocalError:
                        pass
                else:
                    for val in reagent.getchildren():
                        if "identifier" in val.tag:
                            reagent_dict[reagent_type] = val.attrib["value"]
                            break
            if "solvent" not in reagent_dict:
                reagent_dict["solvent"] = ""
            if "catalyst" not in reagent_dict:
                reagent_dict["catalyst"] = ""

            return reagent_dict


def make_df(root):
    """Creates a DataFrame from a list of Reaction objects.

    Args:
        root (xml.etree.ElementTree.Element): The root XML element containing all reactions.

    Returns:
    -------
        pd.DataFrame: A DataFrame containing details of reactions extracted from XML.
    """
    all_rxns = [Reaction(tree) for tree in root]
    reactants = [rx.reactants for rx in all_rxns]
    products = [rx.products for rx in all_rxns]
    reactions = [rx.reaction for rx in all_rxns]
    complete_reactions = [rx.complete_reaction for rx in all_rxns]
    solvents = [rx.solvents for rx in all_rxns]
    catalysts = [rx.catalysts for rx in all_rxns]
    refs = [rx.reference for rx in all_rxns]
    titles = [rx.title for rx in all_rxns]
    paragraphs = [rx.paragraph for rx in all_rxns]
    return pd.DataFrame(
        {
            "REACTANTS": reactants,
            "PRODUCTS": products,
            "RXN": reactions,
            "COMPLETE_RXNS": complete_reactions,
            "SOLVENTS": solvents,
            "CATALYSTS": catalysts,
            "REF": refs,
            "TITLE": titles,
            "PARAGRAPH": paragraphs,
        }
    )


def read_year(year: int):
    """Reads XML files from a specified year, processes the reactions, and saves them to a CSV file.

    Args:
        year (int): The year from which to read the XML files.

    """
    folder = f"grants/{year}"
    all_dfs = []
    for file in os.listdir(folder):
        tree = ET.parse(f"{folder}/{file}")
        root = tree.getroot()
        all_dfs.append(make_df(root))
    df = pd.concat(all_dfs)
    df.to_csv(f"{year}_USPTO.csv")


if __name__ == "__main__":
    read_year(1976)
