import numpy as np
import pandas as pd
from scipy import stats
from rdkit.Chem import RDKFingerprint
from rdkit.Chem import AllChem
from rdkit.Chem import MACCSkeys
from rdkit.Chem import DataStructs
from mordred import Calculator, descriptors
from drug_learning.two_dimensions.Input import base_class as bc
from drug_learning.two_dimensions.Errors import errors as er

class MorganFP(bc.Fingerprint):

    fp_name = "_MorganFP"

    def transform(self):
        super().transform()
        fts = []
        self.mol_names = []
        for mol in self.structures:
            try:
                fp = AllChem.GetMorganFingerprintAsBitVect(mol,2,nBits=2048)
            except:
                continue
            arr = np.zeros((0,), dtype=np.int8)
            DataStructs.ConvertToNumpyArray(fp,arr)
            fts.append(arr)
            self.features = np.array(fts)
            self.mol_names.append(mol.GetProp("_Name"))
        self.columns = [str(i) for i in list(range(self.features.shape[1]))]
        return self.features


class MACCS_FP(bc.Fingerprint):

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
        self.columns = [str(i) for i in list(range(self.features.shape[1]))]
        return self.features


class RDkitFP(bc.Fingerprint):

    fp_name = "_RDkitFP"

    def transform(self):
        super().transform()
        fts = []
        self.mol_names = []
        for mol in self.structures:
            fp = RDKFingerprint(mol)
            arr = np.zeros((0,), dtype=np.int8)
            DataStructs.ConvertToNumpyArray(fp,arr)
            fts.append(arr)
            self.features = np.array(fts)
            self.mol_names.append(mol.GetProp("_Name"))
        self.columns = [str(i) for i in list(range(self.features.shape[1]))]
        return self.features

class UnfoldedRDkitFP(bc.Fingerprint):

    fp_name = "_UnfoldedRDkitFP"

    def __init__(self, voc):
        super().__init__()
        if not voc:
            raise er.NotVocabularyUnfolded("Vocabulary (--voc) must be pass to use unfolded rdkit fingerprints")
        filename, file_extension = os.path.splitext(voc)
        if file_extension != ".npy":
            raise er.IncorrectFormat("Vocabulary must be an npy file (numpy.save)")
        self.voc = np.load(voc)

    def transform(self):
        super().transform()
        fts = []
        self.mol_names = []
        for mol in self.structures:
            fingerprint = []
            fp = AllChem.UnfoldedRDKFingerprintCountBased(mol)
            fpDict = fp.GetNonzeroElements()
            for fragment in self.voc:
                if fragment in fpDict:
                    fingerprint.append(fpDict[fragment])
                else:
                    fingerprint.append(0)
            fts.append(np.array(fingerprint))
            self.mol_names.append(mol.GetProp("_Name"))
        self.features = np.array(fts)
        self.columns = [str(i) for i in list(range(self.features.shape[1]))]
        return self.features


class MordredFP(bc.Fingerprint):

    fp_name = "_MordredFP"

    def transform(self):
        print("Transforming to mordred fp...")
        super().transform()
        self.mol_names = []
        calc = Calculator(descriptors, ignore_3D=True)
        self.df = calc.pandas(self.structures, nproc=1, quiet=True)
        self.columns = self.df.columns
        self.features = self.df.values
        self.mol_names = [mol.GetProp("_Name") for mol in self.structures]
        return self.features

    def clean(self):
        self.df = self.df.apply(pd.to_numeric, errors="coerce")
        self.df = self.df.astype("float64")
        self.df[self.df > 1e3] = None
        self.df = self.df.dropna(axis=1)
        mask = ~np.any(np.abs(stats.zscore(self.df)) < 2, axis=0)
        self.df = self.df.drop(columns=self.df.columns[mask])

        self.columns = self.df.columns
        self.features = self.df.values
        return self.features
