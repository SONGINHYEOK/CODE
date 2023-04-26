import pandas as pd

id_list = ["1344096268", "1344118153"]
sm_list = ["O=C(NCc1cn2cc(C[NH2+]CC3CCCCC3)ccc2n1)c1cc(=O)n2ccccc2n1", "CC1(C)CC[NH+](Cc2ccc(C(=O)NC[C@]3(O)CCCN(c4cc(NCc5ccccc5)ncn4)C3)c(O)c2)CC1"]


df =pd.DataFrame()

df['id'] = id_list
df['smiles'] = sm_list


df.to_csv('./116_mol.csv', index=False)