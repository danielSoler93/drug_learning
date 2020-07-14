import numpy as np
from rdkit.Chem import RDKFingerprint
from rdkit.Chem import AllChem
from rdkit.Chem import MACCSkeys
from rdkit.Chem import DataStructs
from mordred import Calculator, descriptors
from drug_learning.two_dimensions.Input import base_class as bc

class MorganFP(bc.Fingerprint):

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
        return self.features


class MordredFP(bc.Fingerprint):

    fp_name = "_MordredFP"

    def transform(self):
        super().transform()
        self.mol_names = []
        calc = Calculator(descriptors, ignore_3D=True)
        self.features = calc.pandas(self.structures)
        self.mol_names = [mol.GetProp("_Name") for mol in self.structures]
        return self.features
