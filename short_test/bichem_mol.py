import os

from torchdrug import data, utils
from torchdrug.core import Registry as R
from torchdrug.utils import doc


@R.register("datasets.PhotoTox")
@doc.copy_args(data.MoleculeDataset.load_csv, ignore=("smiles_field", "target_fields"))
class BIChem(data.MoleculeDataset):
    """
    
    bichem molecie.
    Statistics:
        - #Molecule: 81
        - #Regression task: 2
    Parameters:
        path (str): path to store the dataset
        verbose (int, optional): output verbose level
        **kwargs
    
    """

    csv_file = "./molecule-datasets/phototox_pih.csv"
    target_fields = ["logP"]

    def __init__(self, path, verbose=1, **kwargs):
        self.load_csv(self.csv_file, smiles_field="smiles", target_fields=self.target_fields,
                      verbose=verbose, **kwargs)

    
