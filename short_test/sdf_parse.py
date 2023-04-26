from rdkit import Chem, RDConfig
from rdkit.Chem import MolStandardize, Draw, AllChem, rdMolAlign

sample = "./sample.sdf"

def sdf_parse(sdf_file):
    suppl = Chem.SDMolSupplier(sdf_file)
    for mol in suppl:
        smiles = Chem.MolToSmiles(mol)
        mol_img = Chem.MolFromSmiles(smiles)
        id = mol.GetProp('_Name')
        mw = mol.GetProp('Molecular Weight')
        tp = mol.GetProp('TPSA')
        logp = mol.GetProp('LogP')
        ds = mol.GetProp('Density')
        return print(id, mol_img, smiles, mw, tp, logp, ds)
    
sdf_parse(sample)
  
    
        
        