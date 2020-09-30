import numpy as np
import os
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
            fp = AllChem.GetMorganFingerprintAsBitVect(mol,2,nBits=2048)
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
        filename, file_extension = os.path.splitext(voc)
        if file_extension != ".npy":
            er.IncorrectFormat("Vocabulary must be an npy file (numpy.save)")
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
        super().transform()
        self.mol_names = []
        calc = Calculator(descriptors, ignore_3D=True)
        mols = []
        for mol in self.structures:
            mols.append(mol)
            self.mol_names.append(mol.GetProp("_Name"))
        df = calc.pandas(mols)
        self.columns = df.columns
        self.features = df.values
        return self.features
