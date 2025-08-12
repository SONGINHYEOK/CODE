from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, Crippen
from .models import Molecule

class MoleculeService:
    @staticmethod
    def calculate_properties(smiles):
        """Calculate molecular properties using RDKit"""
        try:
            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                return None
            
            return {
                'canonical_smiles': Chem.MolToSmiles(mol, canonical=True),
                'inchi': Chem.MolToInchi(mol),
                'inchi_key': Chem.MolToInchiKey(mol),
                'molecular_weight': Descriptors.MolWt(mol),
                'logp': Descriptors.MolLogP(mol),
                'hbd': Descriptors.NumHDonors(mol),
                'hba': Descriptors.NumHAcceptors(mol),
            }
        except Exception as e:
            print(f"Error calculating properties: {e}")
            return None
    
    @staticmethod
    def get_or_create_molecule(smiles):
        """Get or create molecule with calculated properties"""
        properties = MoleculeService.calculate_properties(smiles)
        if not properties:
            return None
        
        molecule, created = Molecule.objects.get_or_create(
            inchi_key=properties['inchi_key'],
            defaults={
                'smiles': smiles,
                **properties
            }
        )
        
        return molecule