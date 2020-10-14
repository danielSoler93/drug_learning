import os
import pandas as pd
from rdkit import Chem
from drug_learning.two_dimensions.Errors import errors as er
from drug_learning.two_dimensions.Output import output as ot


class Mol2MolSupplier():
    def __init__(self, file_input, sanitize=True):
        self.file_input = file_input
        self.sanitize = sanitize

    def __iter__(self):
        with open(self.file_input, 'r') as f:
            line = f.readline()
            while not f.tell() == os.fstat(f.fileno()).st_size:
                if line.startswith("@<TRIPOS>MOLECULE"):
                    mol = []
                    mol.append(line)
                    line = f.readline()
                    while not line.startswith("@<TRIPOS>MOLECULE"):
                        mol.append(line)
                        line = f.readline()
                        if f.tell() == os.fstat(f.fileno()).st_size:
                            mol.append(line)
                            break
                    mol[-1] = mol[-1].rstrip() # removes blank line at file end
                    block = ",".join(mol).replace(',', '')
                    mol = Chem.MolFromMol2Block(block, sanitize=self.sanitize)
                    if mol is not None:
                        yield mol


class Fingerprint(ot.Saver):

    fp_name = ""

    def __init__(self):
        self.filename = None
        self.structures = None
        self.features = None
        self.fitted = False

    def fit(self, input_file):
        (self.filename, ext) = os.path.splitext(input_file)
        if ext == ".mol2":
            self.structures = Mol2MolSupplier(input_file)
        else:
            # assume for the moment that any file that is not mol2 will be an
            # sdf, might want stricter parsing if more formats are to be
            # supported
            self.structures = (mol for mol in Chem.SDMolSupplier(input_file) if mol is not None) 
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
            raise er.NotOutputFormat("Must specify an output format")
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
