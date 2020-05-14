from rdkit import Chem
import drug_learning.two_dimension.Graph.atoms as at
import drug_learning.two_dimension.Graph.bonds as bd

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

    def build_graph(self):
        atomic_features = self.build_atomic_features()
        bond_features = self.build_bond_features()
        


    def build_atomic_features(self):
        all_atoms = []
        for atom in self.mol.GetAtoms(): #Order??
            Z = atom.GetAtomicNum()
            neighbours = atom.GetDegree()
            formal_charge= atom.GetFormalCharge()
            atom = at.Atom(Z, neighbours, formal_charge)
            all_atoms.append(atom)
        return all_atoms

    def build_bond_features(self):
        all_bonds = []
        for atom in self.mol.GetAtoms():  #Order??
            for bond in atom.GetBonds():
                connects = [bond.GetBeginAtom(), bond.GetEndAtom()]
                aromatic = bond.GetIsAromatic()
                in_ring = 1 if bond.IsInRing() else 0
                conjugated = 1 if bond.GetIsConjugated() else 0
                order = bond.GetBondType()
                bond = bd.Bond(order, aromatic, conjugated, in_ring, connects)
                all_bonds.append(bond)
        return all_bonds
