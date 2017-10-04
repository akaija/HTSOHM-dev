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
from htsohm.pseudo_material import PseudoMaterial
from htsohm.db import session, Material

def random_number_density(number_density_limits, lattice_constants):
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
    v = lattice_constants["a"] * lattice_constants["b"] * lattice_constants["c"]
    min_atoms = int(min_ND * v)
    max_atoms = int(max_ND * v)
    if min_atoms < 2:
        min_atoms = int(2)
    atoms = randrange(min_atoms, max_atoms + 1, 1)
    return atoms

def generate_pseudo_material(run_id, number_of_atomtypes):
    """Write .def and .cif files for a randomly-generated porous material.

    Args:
        run_id (str): identification string for run.
        number_of_atomtypes (int): number of different chemical species used to
            populate the unit cell.

    Returns:
        material (sqlalchemy.orm.query.Query): database row for storing 
            simulation data specific to the material. See
            `htsohm/db/material.py` for more information.

    Atom-sites and unit-cell dimensions are stored in a .cif-file within RASPA's
    structures library: `$(raspa-dir)/structures/cif`. Force field parameters
    are stored within that material's directory in the RASPA library: 
    `$(raspa-dir)/forcefield/(run_id)-(uuid)`. Partial charges, atomic radii,
    and more define each chemical species in `pseudo_atoms.def`. Lennard-Jones
    type interactions (sigma, epsilon-values) are defined in 
    `force_field_mixing_rules.def`. These interactions can be overwritten in 
    `force_field.def`, but by default no interactions are overwritten in this file.
 
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
    pseudo_material = PseudoMaterial(material.uuid)
    
    pseudo_material.run_id = run_id

    pseudo_material.atom_types = []
    for chemical_id in range(number_of_atomtypes):
        pseudo_material.atom_types.append({
            "chemical-id" : "A_%s" % chemical_id,
            "charge"      : 0.,    # See NOTE above.
            "epsilon"     : uniform(*epsilon_limits),
            "sigma"       : uniform(*sigma_limits)
        })

    pseudo_material.lattice_constants = {}
    for i in ['a', 'b', 'c']:
        pseudo_material.lattice_constants[i] = uniform(*lattice_limits)

    pseudo_material.lattice_angles = {}
    for i in ['alpha', 'beta', 'gamma']:
        pseudo_material.lattice_angles[i] = 90.

    pseudo_material.number_of_atoms   = random_number_density(
        number_density_limits, pseudo_material.lattice_constants)

    pseudo_material.atom_sites = []
    for atom in range(pseudo_material.number_of_atoms):
        atom_site = {"chemical-id" : choice(pseudo_material.atom_types)["chemical-id"]}
        for i in ['x-frac', 'y-frac', 'z-frac']:
            atom_site[i] = random()
        pseudo_material.atom_sites.append(atom_site)

    return material, pseudo_material

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

def mutate_pseudo_material(parent_material, parent_pseudo_material, mutation_strength, generation):
    """Modifies a "parent" material's definition files by perturbing each
    parameter by some factor, dictated by the `mutation_strength`.
    
    Args:
        run_id (str): identification string for run.
        parent_id (str): uuid, identifying parent material in database.
        generation (int): iteration count for overall bin-mutate-simulate routine.
        mutation_strength (float): perturbation factor [0, 1].

    Returns:
        new_material (sqlalchemy.orm.query.Query): database row for storing 
            simulation data specific to the material. See
            `htsohm/db/material.py` for more information.

    Todo:
        * Add methods for assigning and mutating charges.

    """
    
    ########################################################################
    # load boundaries from config-file
    lattice_limits          = config["lattice_constant_limits"]
    number_density_limits   = config["number_density_limits"]

    ########################################################################
    # create material object
    child_material = Material(parent_material.run_id)
    child_material.parent_id = parent_material.id
    child_material.generation = generation
    child_pseudo_material = PseudoMaterial(child_material.uuid)

    child_pseudo_material.run_id = parent_pseudo_material.run_id

    ########################################################################
    # perturb LJ-parameters
    child_pseudo_material.atom_types = parent_pseudo_material.atom_types.copy()
    for atom_type in child_pseudo_material.atom_types:    
        for x in ['epsilon', 'sigma']:
            old_x = atom_type[x]
            random_x = uniform(*config["{0}_limits".format(x)])
            atom_type[x] += mutation_strength * (random_x - old_x)

    ########################################################################
    # calculate new lattice constants
    child_pseudo_material.lattice_constants = parent_pseudo_material.\
            lattice_constants.copy()
    for i in ['a', 'b', 'c']:
        old_x = parent_pseudo_material.lattice_constants[i]
        random_x = uniform(*lattice_limits)
        child_pseudo_material.lattice_constants[i] += mutation_strength * (
                random_x - old_x)

    child_pseudo_material.lattice_angles = parent_pseudo_material.lattice_angles.copy()

    ########################################################################
    #perturb number density, calculate number of atoms
    child_ND = parent_pseudo_material.number_density()
    random_ND = uniform(*number_density_limits)
    child_ND += mutation_strength * (random_ND - child_ND)
    child_number_of_atoms = (
            int(child_ND * child_pseudo_material.volume()))

    ########################################################################
    # remove excess atom-sites, if any
    child_pseudo_material.atom_sites = np.random.choice(
            parent_pseudo_material.atom_sites,
            min(child_number_of_atoms,
                len(parent_pseudo_material.atom_sites)),
            replace=False).tolist()
    
    ########################################################################
    # perturb atom-site positions
    for atom_site in child_pseudo_material.atom_sites:
        for i in ['x-frac', 'y-frac', 'z-frac']:
            atom_site[i] = random_position(
                    atom_site[i], random(), mutation_strength)

    ########################################################################
    # add atom-sites, if needed
    if child_number_of_atoms > len(child_pseudo_material.atom_sites):
        for new_sites in range(child_number_of_atoms -
                len(child_pseudo_material.atom_sites)):
            new_atom_site = {'chemical-id' : choice(
                child_pseudo_material.atom_types)['chemical-id']}
            for i in ['x-frac', 'y-frac', 'z-frac']:
                new_atom_site[i] = random()
            child_pseudo_material.atom_sites.append(new_atom_site)

    return child_material, child_pseudo_material

def mutate_pseudo_material2(parent_material, parent_pseudo_material, mutation_strength, generation, config):
    """Modifies a "parent" material's definition files by perturbing each
    parameter by some factor, dictated by the `mutation_strength`.
    
    Args:
        run_id (str): identification string for run.
        parent_id (str): uuid, identifying parent material in database.
        generation (int): iteration count for overall bin-mutate-simulate routine.
        mutation_strength (float): perturbation factor [0, 1].

    Returns:
        new_material (sqlalchemy.orm.query.Query): database row for storing 
            simulation data specific to the material. See
            `htsohm/db/material.py` for more information.

    Todo:
        * Add methods for assigning and mutating charges.

    """
    
    ########################################################################
    # load boundaries from config-file
    lattice_limits          = config["lattice_constant_limits"]
    number_density_limits   = config["number_density_limits"]

    ########################################################################
    # create material object
    child_material = Material(parent_material.run_id)
    child_material.parent_id = parent_material.id
    child_material.generation = generation
    child_pseudo_material = PseudoMaterial(child_material.uuid)

    child_pseudo_material.run_id = parent_pseudo_material.run_id

    ########################################################################
    # perturb LJ-parameters
    child_pseudo_material.atom_types = parent_pseudo_material.atom_types.copy()
    for atom_type in child_pseudo_material.atom_types:    
        for x in ['epsilon', 'sigma']:
            old_x = atom_type[x]
            random_x = uniform(*config["{0}_limits".format(x)])
            atom_type[x] += mutation_strength * (random_x - old_x)

    ########################################################################
    # calculate new lattice constants
    child_pseudo_material.lattice_constants = parent_pseudo_material.\
            lattice_constants.copy()
    for i in ['a', 'b', 'c']:
        old_x = parent_pseudo_material.lattice_constants[i]
        random_x = uniform(*lattice_limits)
        child_pseudo_material.lattice_constants[i] += mutation_strength * (
                random_x - old_x)

    child_pseudo_material.lattice_angles = parent_pseudo_material.lattice_angles.copy()

    ########################################################################
    #perturb number density, calculate number of atoms
    child_ND = parent_pseudo_material.number_density()
    random_ND = uniform(*number_density_limits)
    child_ND += mutation_strength * (random_ND - child_ND)

    print('\nNUMBER DENSITY : {}\n'.format(child_ND))

    child_number_of_atoms = (
            int(child_ND * child_pseudo_material.volume()))

    ########################################################################
    # remove excess atom-sites, if any
    child_pseudo_material.atom_sites = np.random.choice(
            parent_pseudo_material.atom_sites,
            min(child_number_of_atoms,
                len(parent_pseudo_material.atom_sites)),
            replace=False).tolist()
    
##    xfrac = [e['x-frac'] for e in parent_pseudo_material.atom_sites]
##    yfrac = [e['y-frac'] for e in parent_pseudo_material.atom_sites]
##    zfrac = [e['z-frac'] for e in parent_pseudo_material.atom_sites]
##    max_frac = 1.5 * max(max(xfrac), max(yfrac), max(zfrac))
##    min_frac = 0.5 * min(min(xfrac), min(yfrac), min(zfrac))

    ########################################################################
    # perturb atom-site positions
    for atom_site in child_pseudo_material.atom_sites:
        for i in ['x-frac', 'y-frac', 'z-frac']:
            random_frac = uniform(-1, 1)
            atom_site[i] += mutation_strength * (random_frac - atom_site[i])
#            atom_site[i] = random_position(
#                    atom_site[i], random(), mutation_strength)

    ########################################################################
    # add atom-sites, if needed
    if child_number_of_atoms > len(child_pseudo_material.atom_sites):
        for new_sites in range(child_number_of_atoms -
                len(child_pseudo_material.atom_sites)):
            new_atom_site = {'chemical-id' : choice(
                child_pseudo_material.atom_types)['chemical-id']}
            for i in ['x-frac', 'y-frac', 'z-frac']:
                new_atom_site[i] = uniform(-1, 1)
            child_pseudo_material.atom_sites.append(new_atom_site)

    return child_material, child_pseudo_material

def write_cif_file(material, simulation_path):
    """Writes .cif file for structural information.

    Args:
        file_name (str): path to .cif-file, for example:
            `$(raspa-dir)/structures/cif/(run_id)-(uuid).cif`.
        lattice_constants (dict): crystal lattice constants, for example:
            {"a" : (float),
             "b" : (float),
             "c" : (float)}
        atom_sites (list): dictionaries for each atom-site describing position
            and chemical species, for example:
            {"chemical-id" : chemical_id,
             "x-frac"      : x,
             "y-frac"      : y,
             "z-frac"      : z}

    Writes .cif file in RASPA's library:
        `$(raspa-dir)/structures/cif/(run_id)-(uuid).cif`

    """
    file_name = os.path.join(simulation_path, '%s.cif' % material.uuid)
    with open(file_name, "w") as cif_file:
        cif_file.write(
            "\nloop_\n" +
            "_symmetry_equiv_pos_as_xyz\n" +
            "  x,y,z\n"
        )
        for i in ['a', 'b', 'c']:
            cif_file.write(
            "_cell_length_{}    {}\n".format(i, 
                round(material.lattice_constants[i], 4))
        )
        print(material.lattice_angles)
        for i in ['alpha', 'beta', 'gamma']:
            cif_file.write(
            "_cell_angle_{}     {}\n".format(i,
                round(material.lattice_angles[i], 4))
        )
        cif_file.write(
            "loop_\n" +
            "_atom_site_label\n" +
            "_atom_site_type_symbol\n" +
            "_atom_site_fract_x\n" +
            "_atom_site_fract_y\n" +
            "_atom_site_fract_z\n"
        )
        for atom_site in material.atom_sites:
            cif_file.write(
            "{0:5} C {1:4f} {2:4f} {3:4f}\n".format(
                atom_site["chemical-id"],
                round(atom_site["x-frac"], 4),
                round(atom_site["y-frac"], 4),
                round(atom_site["z-frac"], 4)
            ))

def write_mixing_rules(material, simulation_path):
    """Writes .def file for forcefield information.

    Args:
        file_name (str): path to mixing-rules file, for example:
            `$(raspa-dir)/forcefield/(run_id)-(uuid)/force_field_mixing_rules.def`
        atom_types (list): dictionaries for each atom-type describing LJ-type
            interactions, for example:
            {"chemical-id" : chemical_ids[i],
             "charge"      : 0.,
             "epsilon"     : epsilon,
             "sigma"       : sigma}

    Writes file within RASPA's library:
        `$(raspa-dir)/forcefield/(run_id)-(uuid)/force_field_mixing_rules.def`

    """
    adsorbate_LJ_atoms = [
            ['N_n2',    36.0,       3.31],
            ['C_co2',   27.0,       2.80],
            ['O_co2',   79.0,       3.05],
            ['CH4_sp3', 158.5,      3.72],
            ['He',      10.9,       2.64],
            ['H_com',   36.7,       2.958],
            ['Kr',      167.06,     3.924],
            ['Xe',      110.704,    3.690]
    ]
 
    adsorbate_none_atoms = ['N_com', 'H_h2']

    file_name = os.path.join(simulation_path, 'force_field_mixing_rules.def')
    with open(file_name, "w") as mixing_rules_file:
        mixing_rules_file.write(
            "# general rule for shifted vs truncated\n" +
            "shifted\n" +
            "# general rule tailcorrections\n" +
            "no\n" +
            "# number of defined interactions\n" +
            "{}\n".format(len(material.atom_types) + 10) +
            "# type interaction, parameters.    " +
            "IMPORTANT: define shortest matches first, so" +
            " that more specific ones overwrites these\n"
        )
        for atom_type in material.atom_types:
            mixing_rules_file.write(
                "{0:12} lennard-jones {1:8f} {2:8f}\n".format(
                    atom_type["chemical-id"],
                    round(atom_type["epsilon"], 4),
                    round(atom_type["sigma"], 4)
                )
            )
        for at in adsorbate_LJ_atoms:
            mixing_rules_file.write(
                "{0:12} lennard-jones {1:8f} {2:8f}\n".format(at[0], at[1], at[2])
            )
        for at in adsorbate_none_atoms:
            mixing_rules_file.write(
                "{0:12} none\n".format(at)
            )
        mixing_rules_file.write(
            "# general mixing rule for Lennard-Jones\n" +
            "Lorentz-Berthelot")

def write_pseudo_atoms(material, simulation_path):
    """Writes .def file for chemical information.

    Args:
        file_name (str): path to pseudo atoms definitions file, for example:
            `$(raspa-dir)/forcefield/(run_id)-(uuid)/pseudo_atoms.def`
        atom_types (list): dictionaries for each chemical species, including
            an identification string and charge, for example:
            interactions, for example:
            {"chemical-id" : chemical_ids[i],
             "charge"      : 0.,
             "epsilon"     : epsilon,
             "sigma"       : sigma}

    Returns:
        Writes file within RASPA's library,
        `$(raspa-dir)/forcefield/(run_id)-(uuid)/pseudo_atoms.def`

        NOTE: ALL CHARGES ARE 0. IN THIS VERSION.

    """
    temporary_charge = 0.

    file_name = os.path.join(simulation_path, 'pseudo_atoms.def')
    with open(file_name, "w") as pseudo_atoms_file:
        pseudo_atoms_file.write(
            "#number of pseudo atoms\n" +
            "%s\n" % (len(material.atom_types) + 10) +
            "#type  print   as  chem    oxidation   mass    charge  polarization    B-factor    radii   " +
                 "connectivity  anisotropic anisotrop-type  tinker-type\n")
        for atom_type in material.atom_types:
            pseudo_atoms_file.write(
                "{0:7}  yes  C   C   0   12.0       {0:8}  0.0  1.0  1.0    0  0  absolute  0\n".format(
                    atom_type['chemical-id'], 0.0
                )
            )
        pseudo_atoms_file.write(
            "N_n2     yes  N   N   0   14.00674   -0.4048   0.0  1.0  0.7    0  0  relative  0\n" +
            "N_com    no   N   -   0    0.0        0.8096   0.0  1.0  0.7    0  0  relative  0\n" +
            "C_co2    yes  C   C   0   12.0        0.70     0.0  1.0  0.720  0  0  relative  0\n" +
            "O_co2    yes  O   O   0   15.9994    -0.35     0.0  1.0  0.68   0  0  relative  0\n" +
            "CH4_sp3  yes  C   C   0   16.04246    0.0      0.0  1.0  1.00   0  0  relative  0\n" +
            "He       yes  He  He  0    4.002602   0.0      0.0  1.0  1.0    0  0  relative  0\n" +
            "H_h2     yes  H   H   0    1.00794    0.468    0.0  1.0  0.7    0  0  relative  0\n" +
            "H_com    no   H   H   0    0.0        0.936    0.0  1.0  0.7    0  0  relative  0\n" +
            "Xe       yes  Xe  Xe  0  131.293      0.0      0.0  1.0  2.459  0  0  relative  0\n" +
            "Kr       yes  Kr  Kr  0   83.798      0.0      0.0  1.0  2.27   0  0  relative  0\n"
        )

def write_force_field(simulation_path):
    """Writes .def file to overwrite LJ-type interactions.

    Args:
        file_name (str): path to .def-file, for example:
            `$(raspa-dir)/forcefield/(run_id)-(uuid)/force_field.def`

    Writes file within RASPA's library:
        `$(raspa-dir)/forcefield/(run_id)-(uuid)/force_field.def`

    NOTE: NO INTERACTIONS ARE OVERWRITTEN BY DEFAULT.

    """
    file_name = os.path.join(simulation_path, 'force_field.def')
    with open(file_name, "w") as force_field_file:
        force_field_file.write(
            "# rules to overwrite\n" +
            "0\n" +
            "# number of defined interactions\n" +
            "0\n" +
            "# mixing rules to overwrite\n" +
            "0")
