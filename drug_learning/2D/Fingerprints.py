import os
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs

class MorganFP():
    def __init__(self):
        self.structures = None
        self.features = None

    def fit(self, input_sdf):
        if self.structures is None:
            self.structures = Chem.SDMolSupplier(input_sdf)
        return self.structures

    def transform(self, to_csv = False):
        fts = []
        fps = []
        if self.features is None:
            for mol in self.structures:
                fp = AllChem.GetMorganFingerprintAsBitVect(mol,2,nBits=1024)
                fps.append(fp)
                arr = np.zeros((0,), dtype=np.int8)
                DataStructs.ConvertToNumpyArray(fp,arr)
                fts.append(arr)
                self.features = np.array(fts)
        return self.features


if __name__ == "__main__":

    morgan_fps = MorganFP()
    morgan_fps.fit("../datasets/ligands.sdf")
    morgan_fps.structures
    morgan_fps.transform()
    morgan_fps.features
