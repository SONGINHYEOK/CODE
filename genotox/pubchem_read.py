import pandas as pd 

df = pd.read_csv("/Users/song-inhyeog/Downloads/PubChem_bioassay_text_Genotoxicity.csv")

col_list  = df.columns

print(df["aidsrcname"])

