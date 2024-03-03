import xml.etree.ElementTree as ET
import pandas as pd
import numpy as np
from tqdm import tqdm
import os


class Reaction:
    def __init__(self, tree):
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
        return self.tree[0].getchildren()[0].text
    
    def get_rxn_smiles(self):
        return self.tree[1].text
    
    def get_titles(self):
        return self.tree[0][1].text
    
    def get_paragraphs(self):
        return self.tree[0][-1].text
    
    def get_reaction(self):
        comps = self.complete_reaction.split(">")
        return comps[0], comps[-1], f"{comps[0]}>>{comps[-1]}"
    
    def get_reagents(self):
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
                        reagent_dict[reagent_type] = f"{reagent_dict[reagent_type]}.{reag_smile}"
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
    return pd.DataFrame({"REACTANTS": reactants, "PRODUCTS": products, 
                         "RXN": reactions, "COMPLETE_RXNS": complete_reactions, 
                         "SOLVENTS": solvents, "CATALYSTS": catalysts, "REF": refs, 
                         "TITLE": titles, "PARAGRAPH": paragraphs})
    

def read_year(year: int):
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
