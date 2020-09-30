import os
import pandas as pd
from rdkit import Chem
from drug_learning.two_dimensions.Errors import errors as er
from drug_learning.two_dimensions.Output import output as ot

class Fingerprint(ot.Saver):

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
            raise er.NotFittedException("Must fit the model before transform")

    def clean(self):
        return self.features

    def save(self, to_csv=False, to_parquet=False, to_feather=False, to_hdf=False,
            to_pickle=False):
        if not to_csv and not to_parquet and not to_feather and not to_hdf and not  to_pickle:
            raise er.NotFormat("Must specify an output format")
        if not self.mol_names:
            raise er.NotTransformException("Must transform the input molecules before save")
        self.df = pd.DataFrame(self.features, index=self.mol_names, columns=self.columns)
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
