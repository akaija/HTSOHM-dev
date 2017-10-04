# standard library imports
import sys
import os
from random import choice, random, randrange, uniform
import shutil
from uuid import uuid4

# related third party imports
import numpy as np
import yaml

# local application/library specific imports
import htsohm
from htsohm import config
from htsohm.db import session, Material, Structure, LennardJones, AtomSites

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

def generate_material(run_id, number_of_atomtypes):
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
    material.generation = 0

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

def closest_distance(x, y):
    """Finds closest distance between two points across periodic boundaries.

    Args:
        x (float): value in range [0,1].
        y (float): value in range [0,1].

    Returns:
        D (float): closest distance measured between `x` and `y`, with
            periodic boundaries at 0 and 1. For example, the disance between
            x = 0.2 and y = 0.8 would be 0.4, and not 0.6.

    """
    a = 1 - y + x
    b = abs(y - x)
    c = 1 - x + y
    return min(a, b, c)

def random_position(x_o, x_r, strength):
    """Produces a point along the closest path between two points.

    Args:
        x_o (float): value in range [0,1].
        x_r (float): value in range [0,1].
        strength (float): refers to mutation strength. Value determines
            fractional distance to `x_r` from `x_o`.

    Returns:
        xfrac (float): a value between `x_o` and `x_r` across the closest
            distance accross periodic boundaries at 0 and 1.

    """
    dx = closest_distance(x_o, x_r)
    if (x_o > x_r
            and (x_o - x_r) > 0.5):
        xfrac = (x_o + strength * dx) % 1.
    if (x_o < x_r
            and (x_r - x_o) > 0.5):
        xfrac = (x_o - strength * dx) % 1.
    if (x_o >= x_r
            and (x_o - x_r) < 0.5):
        xfrac = x_o - strength * dx
    if (x_o < x_r
            and (x_r - x_o) < 0.5):
        xfrac = x_o + strength * dx
    return xfrac

def mutate_material(parent_material, mutation_strength, generation):
    """Modifies a "parent" pseudomaterial structure by perturbing each
    parameter by some factor, dictated by the `mutation_strength`.
    
    Args:
        parent_material : record for parent pseudomaterial in 'materials' table.
        mutation_strength (float): perturbation factor [0, 1].
        generation (int): iteration count for overall bin-mutate-simulate routine.

    Returns:

    Todo:
        * Add methods for assigning and mutating charges.

    """
    
    ########################################################################
    # load boundaries from config-file
    lattice_limits          = config["lattice_constant_limits"]
    number_density_limits   = config["number_density_limits"]

    child_material = Material(parent_material.run_id)
    child_material.generation = generation
    child_material.parent_id = parent_material.id

    print('Parent UUID :\t{}'.format(parent_material.uuid))
    print('Child UUID :\t{}'.format(child_material.uuid))

    # perturb lennard-jones parameters
    for atom_type in parent_material.structure.lennard_jones:
        child_material.structure.lennard_jones.append(LennardJones(
            chemical_id = atom_type.chemical_id,
            sigma = atom_type.sigma + mutation_strength * (uniform(*config['sigma_limits']) - atom_type.sigma),
            epsilon = atom_type.epsilon + mutation_strength * (uniform(*config['epsilon_limits']) - atom_type.epsilon)))

    if config['interactive_mode'] == 'on':
        print('=========================================================================================')
        print('LENNARD JONES PARAMETERS')
        print('=========================================================================================')
        print('  chemical-id\t|  parent sigma\t|  child sigma\t|  parent epsilon\t| child epsilon')
        print('----------------+---------------+---------------+-----------------------+----------------')
        for i in range(len(child_material.structure.lennard_jones)):
            c_chem = child_material.structure.lennard_jones[i]
            p_chem = parent_material.structure.lennard_jones[i]
            print('  {}\t\t|  {}\t|  {}\t|  {}\t\t|  {}'.format(c_chem.chemical_id,
                round(p_chem.sigma, 4), round(c_chem.sigma, 4),
                round(p_chem.epsilon, 4), round(c_chem.epsilon, 4)))

    # perturb lattice constants
    child_material.structure.lattice_constant_a = parent_material.structure.lattice_constant_a \
            + mutation_strength * (uniform(*lattice_limits) - parent_material.structure.lattice_constant_a)
    child_material.structure.lattice_constant_b = parent_material.structure.lattice_constant_b \
            + mutation_strength * (uniform(*lattice_limits) - parent_material.structure.lattice_constant_b)
    child_material.structure.lattice_constant_c = parent_material.structure.lattice_constant_c \
            + mutation_strength * (uniform(*lattice_limits) - parent_material.structure.lattice_constant_c)

    if config['interactive_mode'] == 'on':
        print('==========================================')
        print('LATTICE CONSTANTS')
        print('==========================================')
        print('  direction\t|  parent \t|  child')
        print('----------------+---------------+---------')
        for i in ['a', 'b', 'c']:
            print('  {}\t\t|  {}\t|  {}'.format(i,
                round(getattr(parent_material.structure, 'lattice_constant_{}'.format(i)), 4),
                round(getattr(child_material.structure, 'lattice_constant_{}'.format(i)), 4)))

    # perturb number density/number of atom-sites
    parent_ND = len(parent_material.structure.atom_sites) / parent_material.structure.volume
    child_ND = parent_ND + mutation_strength * (uniform(*number_density_limits) - parent_ND)
    number_of_atoms = int(child_ND * child_material.structure.volume)

    # remove atom-sites, if necessary
    child_material.structure.atom_sites = np.random.choice(
        parent_material.structure.atom_sites,
        min(number_of_atoms, len(parent_material.structure.atom_sites)),
        replace = False).tolist()

    # store original atom-site positions before perturbation in order to compare
    # parent and child values later
    if config['interactive_mode'] == 'on':
        p_x, p_y, p_z = [], [], []
        for atom_site in child_material.structure.atom_sites:
            p_x.append(atom_site.x_frac)
            p_y.append(atom_site.y_frac)
            p_z.append(atom_site.z_frac)

    # perturb atom-site positions
    for atom_site in child_material.structure.atom_sites:
        atom_site.x_frac = random_position(atom_site.x_frac, random(), mutation_strength)
        atom_site.y_frac = random_position(atom_site.y_frac, random(), mutation_strength)
        atom_site.z_frac = random_position(atom_site.z_frac, random(), mutation_strength)

    # add atom-sites, if needed
    if number_of_atoms > len(parent_material.structure.atom_sites):
        for new_sites in range(number_of_atoms - len(parent_material.structure.atom_sites)):
            parent_material.structure.atom_sites.append(AtomSites(
                chemical_id = 'A_{}'.format(choice(
                    range(len(parent_material.structure.lennard_jones)))),
                x_frac = random(), y_frac = random(), z_frac = random()))

    if config['interactive_mode'] == 'on':
        print('===================================================================')
        print('FIRST 10 ATOM-SITES')
        print('===================================================================')
        print('  chemical-id\t|  parent position\t\t|  child position')
        print('----------------+-------------------------------+------------------')
        for i in range(min([10, len(p_x)])):
            c = child_material.structure.atom_sites[i]
            print('  {}\t\t|  {}\t|  {}'.format(c.chemical_id,
                (round(p_x[i], 4), round(p_y[i], 4), round(p_z[i], 4)),
                (round(c.x_frac, 4), round(c.y_frac, 4), round(c.z_frac, 4))))

        print('==========================')
        print('NUMBER DENSITY')
        print('==========================')
        print('  parent\t|  child')
        print('----------------+---------')
        print('  {}\t| {}'.format(
            round(len(parent_material.structure.atom_sites) / parent_material.structure.volume, 6),
            round(len(child_material.structure.atom_sites) / child_material.structure.volume, 6)))

    return child_material
