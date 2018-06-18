from random import choice, random, randrange, uniform

from htsohm import config
from htsohm.db import Material, Structure, LennardJones, AtomSites

def random_number_density(number_density_limits, material_structure): #lattice_constants):
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
    v = material_structure.volume
    min_atoms = int(min_ND * v)
    max_atoms = int(max_ND * v)
    if min_atoms < 2:
        min_atoms = int(2)
    atoms = randrange(min_atoms, max_atoms + 1, 1)
    return atoms

def generate_material(run_id, gen, number_of_atomtypes):
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
    lattice_limits          = config["lattice_constant_limits"]
    number_density_limits   = config["number_density_limits"]
    epsilon_limits          = config["epsilon_limits"]
    sigma_limits            = config["sigma_limits"]
    max_charge              = config["charge_limit"]
    elem_charge             = config["elemental_charge"]

    ########################################################################
    material = Material(run_id)
    material.generation = gen

    structure = material.structure
    structure.lattice_constant_a = uniform(*lattice_limits)
    structure.lattice_constant_b = uniform(*lattice_limits)
    structure.lattice_constant_c = uniform(*lattice_limits)

    for chemical_id in range(number_of_atomtypes):
        structure.lennard_jones.append(LennardJones(
                    chemical_id = 'A_{}'.format(chemical_id),
                    sigma = uniform(*sigma_limits),
                    epsilon = uniform(*epsilon_limits)))

    atom_sites = structure.atom_sites
    for i in range(random_number_density(number_density_limits, structure)):
        atom_sites.append(AtomSites(
            chemical_id = 'A_{}'.format(choice(range(number_of_atomtypes))),
            x_frac = random(), y_frac = random(), z_frac = random()))

    return material


def new_material(run_id, gen):
    return generate_material(run_id, gen, config['number_of_atom_types'])
