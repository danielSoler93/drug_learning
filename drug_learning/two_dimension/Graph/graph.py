from rdkit import Chem
import pandas as pd
import drug_learning.two_dimension.Graph.atoms as at
import drug_learning.two_dimension.Graph.bonds as bd
import drug_learning.two_dimension.Graph.custom_errors as ce

bos={Chem.BondType.SINGLE:1.0,
     Chem.BondType.DOUBLE:2.0,
     Chem.BondType.TRIPLE:3.0,
     Chem.BondType.AROMATIC:1.5,
     Chem.BondType.UNSPECIFIED:0.0}

"""
Try to reproduce methodology fom https://sci-hub.im/10.1021/acs.jcim.6b00601
"""

class GraphBuilder:

    def __init__(self, smiles):
        self.smiles = smiles
        self.mol = Chem.MolFromSmiles(self.smiles)
        self.atomic_features = None
        self.bond_features = None
        self.molecule_features = None

    def build_graph(self):
        if self.mol is None:
            return None
        self.atomic_features = self.build_atomic_features()
        self.bond_features = self.build_bond_features()
        self.molecule_features = self.join_bond_atomic_features()
        return self.molecule_features




    def build_atomic_features(self):
        all_atoms = []
        for atom in self.mol.GetAtoms(): #Order??
            Z = atom.GetAtomicNum()
            neighbours = atom.GetDegree()
            formal_charge= atom.GetFormalCharge()
            index = atom.GetIdx()
            atom = at.Atom(Z, neighbours, formal_charge, index)
            all_atoms.append(atom)
        return all_atoms

    def build_bond_features(self):
        all_bonds = []
        for atom in self.mol.GetAtoms():  #Order??
            for bond in atom.GetBonds():
                connects = [bond.GetBeginAtom().GetIdx(), bond.GetEndAtom().GetIdx()]
                aromatic = 1 if bond.GetIsAromatic() else 0
                in_ring = 1 if bond.IsInRing() else 0
                conjugated = 1 if bond.GetIsConjugated() else 0
                order = bos[bond.GetBondType()]
                bond = bd.Bond(order, aromatic, conjugated, in_ring, connects)
                all_bonds.append(bond)
        return all_bonds

    #Fix to meet paper
    def join_bond_atomic_features(self):

        if not self.atomic_features:
            raise ce.MissAtomicFeatures("Missing atomic features. Call graph.build_atomic_features() first.")

        if not self.bond_features:
            raise ce.MissBondFeatures("Missing bond features. Call graph.build_bond_features() first.")

        molecule_features = [0]*len(self.atomic_features)
        for bond in self.bond_features:
            for atom in bond.connects:
                atomic_features = self.atomic_features[atom].get_vector()
                bond_features = bond.get_vector()
                molecule_features[atom] = atomic_features + bond_features
        return molecule_features