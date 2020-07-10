import os
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs

class NotFittedException(Exception):
    pass

class NotTransformException(Exception):
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
class MorganFP(Saver):
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

    def transform(self, to_csv = False, to_parquet=False):
        if not self.fitted:
            raise NotFittedException("Must fit the model before transform")
        fts = []
        fps = []
        self.mol_names = []
        for mol in self.structures:
            fp = AllChem.GetMorganFingerprintAsBitVect(mol,2,nBits=1024)
            fps.append(fp)
            arr = np.zeros((0,), dtype=np.int8)
            DataStructs.ConvertToNumpyArray(fp,arr)
            fts.append(arr)
            self.features = np.array(fts)
            self.mol_names.append(mol.GetProp("_Name"))

    def save(self, to_csv=False, to_parquet=True, to_feather=False, to_hdf=False,
            to_pickle=False):
        if not self.mol_names:
            raise NotTransformException("Must transform the input molecules before save")
        column_names = [str(i) for i in list(range(self.features.shape[1]))]
        df = pd.DataFrame(self.features, index=self.mol_names, columns=column_names)
        Saver.__init__(self, df)
        if to_csv:
            self.to_csv(self.filename + "_MorganFP.csv")
        if to_parquet:
            self.to_parquet(self.filename + "_MorganFP.gzip")
        if to_feather:
            self.to_feather(self.filename + "_MorganFP.ftr")
        if to_hdf:
            self.to_hdf(self.filename + "_MorganFP.hdf5")
        if to_pickle:
            self.to_pickle(self.filename + "_MorganFP.pkl")
        return self.features


if __name__ == "__main__":

    morgan_fps = MorganFP()
    morgan_fps.fit("../datasets/ligands.sdf")
    morgan_fps.transform()
    morgan_fps.save(to_csv=False, to_parquet=False, to_feather=False, to_hdf=False, to_pickle=True)

