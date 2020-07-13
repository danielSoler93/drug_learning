import os
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import MACCSkeys
from rdkit.Chem import DataStructs
from mordred import Calculator, descriptors

class NotFittedException(Exception):
    pass

class NotFittedException(Exception):
    pass

class Saver():
    def __init__(self, df):
        self.df = df

    def to_csv(self, output):
        self.df.to_csv(output)

    def to_parquet(self, output):
        self.df.to_parquet(output, compression='gzip')

    def to_feather(self, output):
        self.df = self.df.reset_index() #Requiered by pandas. Parquet is a better option
        self.df.to_feather(output)

    def to_hdf(self, output):
        self.df.to_hdf(output, key='df', mode='w')

    def to_pickle(self, output):
        self.df.to_pickle(output)

class Fingerprint(Saver):
    fp_name = ""

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

    def transform(self):
        if not self.fitted:
            raise NotFittedException("Must fit the model before transform")

    def save(self, to_csv=False, to_parquet=True, to_feather=False, to_hdf=False,
            to_pickle=False):
        if not self.mol_names:
            raise NotTransformException("Must transform the input molecules before save")
        column_names = [str(i) for i in list(range(self.features.shape[1]))]
        df = pd.DataFrame(self.features, index=self.mol_names, columns=column_names)
        Saver.__init__(self, df)
        if to_csv:
            self.to_csv(self.filename + self.fp_name + ".csv")
        if to_parquet:
            self.to_parquet(self.filename + self.fp_name + ".gzip")
        if to_feather:
            self.to_feather(self.filename + self.fp_name + ".ftr")
        if to_hdf:
            self.to_hdf(self.filename + self.fp_name + ".hdf5")
        if to_pickle:
            self.to_pickle(self.filename + self.fp_name + ".pkl")
        return self.features


class MorganFP(Fingerprint):

    fp_name = "_MorganFP"

    def transform(self):
        super().transform()
        fts = []
        self.mol_names = []
        for mol in self.structures:
            fp = AllChem.GetMorganFingerprintAsBitVect(mol,2,nBits=1024)
            arr = np.zeros((0,), dtype=np.int8)
            DataStructs.ConvertToNumpyArray(fp,arr)
            fts.append(arr)
            self.features = np.array(fts)
            self.mol_names.append(mol.GetProp("_Name"))
        return self.features


class MACCS_FP(Fingerprint):

    fp_name = "_MACCS_FP"

    def transform(self):
        super().transform()
        fts = []
        self.mol_names = []
        for mol in self.structures:
            fp = MACCSkeys.GenMACCSKeys(mol)
            arr = np.zeros((0,), dtype=np.int8)
            DataStructs.ConvertToNumpyArray(fp,arr)
            fts.append(arr)
            self.features = np.array(fts)
            self.mol_names.append(mol.GetProp("_Name"))
        return self.features


class RDkitFP(Fingerprint):

    fp_name = "_RDkitFP"

    def transform(self):
        super().transform()
        fts = []
        self.mol_names = []
        for mol in self.structures:
            fp = Chem.RDKFingerprint(mol)
            arr = np.zeros((0,), dtype=np.int8)
            DataStructs.ConvertToNumpyArray(fp,arr)
            fts.append(arr)
            self.features = np.array(fts)
            mol_names.append(mol.GetProp("_Name"))
        return self.features


class MordredFP(Fingerprint):

    self.fp_name = "_MordredFP"

    def transform(self):
        super().transform()
        self.mol_names = []
        calc = Calculator(descriptors, ignore_3D=True)
        self.features = calc.pandas(self.structures)
        self.mol_names = [mol.GetProp("_Name") for mol in self.structures]
        return self.features


if __name__ == "__main__":

    morgan_fps = MorganFP()
    morgan_fps.fit("../datasets/ligands.sdf")
    morgan_fps.structures
    morgan_fps.transform()
    morgan_fps.save(to_csv=False, to_parquet=False, to_feather=True, to_hdf=False, to_pickle=True)
