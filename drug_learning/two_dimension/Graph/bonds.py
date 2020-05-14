



class Bond:

    def __init__(self, order, aromatic, conjugated, in_ring, connects):
        self.order = order
        self.aromatic = aromatic
        self.conjugated = conjugated
        self.in_ring = in_ring
        self.connects = connects

    def get_vector(self):
        return [self.order, self.aromatic, self.conjugated, *self.connects]