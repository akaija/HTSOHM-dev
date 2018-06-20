class LatticeConstants:
    def __init__(self, a, b, c):
        self.a            = a
        self.b            = b
        self.c            = c

class AtomSite:
    def __init__(self, chemical_id, x, y, z):
        self.chemical_id  = chemical_id
        self.x            = x
        self.y            = y
        self.z            = z
        self.q            = 0.0

class AtomType:
    def __init__(self, chemical_id, sigma, epsilon):
        self.chemical_id  = chemical_id
        self.sigma        = sigma
        self.epsilon      = epsilon

class Structure:
    """Class for storing pseudomaterial structural data.
    Attributes:
        lattice_constants (LatticeConstants) : crystal lattice parameters.
        atom_types (list: AtomType) : Lennard Jones parameters.
        atom_sites (list: AtomSite) : atom-site positions with partial charges.
    """

    def __init__(self):
        self.lattice_constants = LatticeConstants(None, None, None)
        self.atom_types= []
        self.atom_sites = []

    def __repr__(self):
        return ('{0.__class__.__name__!s}('
                '{0.lattice_constants!r}, '
                '{0.atom_types!r}, '
                '{0.atom_sites!r})').format(self)

    def volume(self):
        return (self.lattice_constants.a *
                self.lattice_constants.b *
                self.lattice_constants.c)

    def n (self):
        return len(self.atom_sites)

    def number_density(self):
        return self.n() / self.volume()
