from math import ceil
from random import randint

from Crypto.Hash import SHA256
from zlib import crc32

from htsohm.db import Material
from htsohm.structure import Structure, LatticeConstants, AtomSite, AtomType

def get_n_digit_seed(n):
    range_start = 10 ** (n - 1)
    range_end = (10 ** n) - 1
    return randint(range_start, range_end)

# random number generator
def bytes_to_float(b):
    return float(crc32(b) & 0xffffffff) / 2**32

def str_to_float(s, encoding="utf-8"):
    return bytes_to_float(s.encode(encoding))

def random_number(n):
    # use SHA256 to get random number from some input value
    h = SHA256.new('{}'.format(n).encode())
    return str_to_float(h.hexdigest())

# uniform selection from random number
def uniform_selection(a, b, r):
    return b - r * (b - a)

def random_number_density(number_density_limits, structure, r): #lattice_constants):
    """Produces random number for atom-sites in a unit cell, constrained by
    some number density limits.

    Args:
        number_density_limits (list): min and max number densities, for example:
            [min_(float), max_(float)]
        lattice_constants (dict): crystal lattice constants, for example:
            {"a" : (float),
             "b" : (float),
             "c" : (float)}
    
    Returns:
        atoms (int): some random number of atom-sites under the imposed limits.

        If the minimum number density results in a unit cell with less than 2
        atom-sites with the given lattice constants, a minimum number density
        of TWO ATOM SITES PER UNIT CELL is imposed.
    
    """
    min_ND = number_density_limits[0]
    max_ND = number_density_limits[1]
    v = structure.volume()
    min_atoms = ceil(min_ND * v)
    max_atoms = int(max_ND * v)
    if min_atoms < 2:
        min_atoms = 2
    return int(uniform_selection(min_atoms, max_atoms, r))

def generate_material(run_id, seed, config):
    """Create records for pseudomaterial simulation and structure data."

    Args:
        run_id (str): identification string for run.
        number_of_atomtypes (int): number of different chemical species used to
            populate the unit cell.

    Returns:
        material (sqlalchemy.orm.query.Query): database row for storing 
            simulation data specific to the material. See
            `htsohm/db/material.py` for more information.

    """
    number_of_atom_types    = config["number_of_atom_types"]
    lattice_limits          = config["lattice_constant_limits"]
    number_density_limits   = config["number_density_limits"]
    epsilon_limits          = config["epsilon_limits"]
    sigma_limits            = config["sigma_limits"]
    max_charge              = config["charge_limit"]

    ########################################################################
    # create database row
    material = Material(run_id, seed)

    # create structure object
    structure = Structure()

    # assign lattice constants
    structure.lattice_constants = LatticeConstants(
            uniform_selection(*lattice_limits, random_number(seed)),
            uniform_selection(*lattice_limits, random_number(seed + 1)),
            uniform_selection(*lattice_limits, random_number(seed + 2)))
    seed += 3

    # store unit cell volume to row
    material.unit_cell_volume = structure.volume()
    
    # assign Lennard-Jones parameters
    for chemical_id in range(number_of_atom_types):
        structure.atom_types.append(
                AtomType("A_{}".format(chemical_id),
                    uniform_selection(*sigma_limits, random_number(seed)),
                    uniform_selection(*epsilon_limits, random_number(seed + 1))))
        seed += 2
    
    # calculate random number of atom-sites
    number_of_atoms = random_number_density(number_density_limits, structure, random_number(seed))
    print("Number of atoms : {}".format(number_of_atoms))
    seed += 1

    # store number density to row
    material.number_density = number_of_atoms / structure.volume()

    # assign atom-site positions and calculate avg. sigma/epsilon values
    sigma_sum, epsilon_sum = 0, 0
    for i in range(number_of_atoms):
        # select chemical species
        chemical_id = int(uniform_selection(0, number_of_atom_types, random_number(seed)))
        seed += 1

        # sum LJ parameters for averaging later
        sigma_sum += structure.atom_types[chemical_id].sigma
        epsilon_sum += structure.atom_types[chemical_id].epsilon

        # select position sufficiently distanced from all other atom-sites
        while True:
            x, y, z = random_number(seed), random_number(seed + 1), random_number(seed + 2)
            seed += 3
            distance_passed = True
            for site in structure.atom_sites:
                distance_squared = (x - site.x) ** 2 + (y - site.y) ** 2 + (z - site.z) ** 2
                if distance_squared <= 10 ** -4:
                    distance_passed = False
                    break
            if distance_passed == True:
                break
        structure.atom_sites.append(AtomSite("A_{}".format(chemical_id), x, y, z))
    material.average_sigma = sigma_sum / structure.n()
    material.average_epsilon = epsilon_sum / structure.n()
    
    # assign atom-site partial charges
    for i in range(structure.n()):
        a0 = max_charge - abs(structure.atom_sites[i].q)
        j = int(uniform_selection(0, structure.n() - 1, random_number(seed)))
        a1 = max_charge - abs(structure.atom_sites[j].q)
        dq = float("{0:.6f}".format(
            uniform_selection(0, min([a0, a1]), random_number(seed + 1))))
        structure.atom_sites[i].q += dq
        structure.atom_sites[j].q -= dq
        seed += 2

    net_charge = 0
    for i in range(structure.n()):
        net_charge += structure.atom_sites[i].q
    print("FRAMEWORK NET CHARGE :\t{}".format(net_charge))

    return material, structure

def new_material(run_id, config):
    # number of digits for seed
    n = 9
    # get seed value (initial value to be incremented later)
    seed = get_n_digit_seed(n)
    print("CREATING MATERIAL WITH SEED : {}".format(seed))
    return generate_material(run_id, seed, config)
