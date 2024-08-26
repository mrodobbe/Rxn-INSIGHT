from rxn_insight.database import *


def test_database_creation():
    db = Database()
    df = pd.DataFrame({"RXN": ["CCOC(=O)c1ccccc1.CNc1ccccc1>>CN(C(=O)c1ccccc1)c1ccccc1",
                               "C#Cc1ccccc1.O=CCCc1ccccc1.C1COCCN1>>C(#CC(CCc1ccccc1)N1CCOCC1)c1ccccc1"]})
    db.create_database_from_df(df=df, reaction_column="RXN")
