import os
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import MACCSkeys
from rdkit.Chem import DataStructs

class NotFittedException(Exception):
    pass

class Fingerprint():
    def __init__(self):
        self.filename = None
        self.structures = None
        self.features = None
        self.fitted = False

    def fit(self, input_sdf):
        (self.filename, ext) = os.path.splitext(input_sdf)
        self.structures = Chem.SDMolSupplier(input_sdf)
        self.fitted = True
        return self.structures

    def transform(self, to_csv = False):
        if not self.fitted:
            raise NotFittedException("Must fit the model before transform")


class MorganFP(Fingerprint):

    def transform(self, to_csv = False):
        super().transform(to_csv = False)
        fts = []
        mol_names = []
        for mol in self.structures:
            fp = AllChem.GetMorganFingerprintAsBitVect(mol,2,nBits=1024)
            arr = np.zeros((0,), dtype=np.int8)
            DataStructs.ConvertToNumpyArray(fp,arr)
            fts.append(arr)
            self.features = np.array(fts)
            mol_names.append(mol.GetProp("_Name"))

        if to_csv:
            df = pd.DataFrame(self.features, index = mol_names)
            df.to_csv(self.filename + "_MorganFP.csv")

        return self.features

class MACCS_FP(Fingerprint):

    def transform(self, to_csv = False):
        super().transform(to_csv = False)
        fts = []
        mol_names = []
        for mol in self.structures:
            fp = MACCSkeys.GenMACCSKeys(mol)
            arr = np.zeros((0,), dtype=np.int8)
            DataStructs.ConvertToNumpyArray(fp,arr)
            fts.append(arr)
            self.features = np.array(fts)
            mol_names.append(mol.GetProp("_Name"))

        if to_csv:
            df = pd.DataFrame(self.features, index = mol_names)
            df.to_csv(self.filename + "_MACCS_FP.csv")

        return self.features


class RDkitFP(Fingerprint):

    def transform(self, to_csv = False):
        super().transform(to_csv = False)
        fts = []
        mol_names = []
        for mol in self.structures:
            fp = Chem.RDKFingerprint(mol)
            arr = np.zeros((0,), dtype=np.int8)
            DataStructs.ConvertToNumpyArray(fp,arr)
            fts.append(arr)
            self.features = np.array(fts)
            mol_names.append(mol.GetProp("_Name"))

        if to_csv:
            df = pd.DataFrame(self.features, index = mol_names)
            df.to_csv(self.filename + "_RDkitFP.csv")

        return self.features


if __name__ == "__main__":

    morgan_fps = MorganFP()
    morgan_fps.fit("../datasets/ligands.sdf")
    morgan_fps.structures
    morgan_fps.transform(to_csv = True)
    morgan_fps.features
