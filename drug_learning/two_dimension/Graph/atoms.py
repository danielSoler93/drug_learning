


class Atom:

    def __init__(self, Z, neighbours, formal_charge, index ):
        self.Z = Z
        self.neighbours = neighbours
        self.formal_charge = formal_charge
        self.index = index

    def get_vector(self):
        return [self.Z, self.neighbours, self.formal_charge]