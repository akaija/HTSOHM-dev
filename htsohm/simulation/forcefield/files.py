# standard library imports
import sys
import os
from random import choice, random, randrange, uniform
import shutil
from uuid import uuid4
from string import Template

# related third party imports
import numpy as np
import yaml

# local application/library specific imports
import htsohm
from htsohm import config
from htsohm.db import session, Material, Structure, LennardJones, AtomSites

def write_cif_file(material, simulation_path):
    """Writes .cif file for structural information.

    Args:


    """
    # Load .cif-file template and replace lattice constant values
    htsohm_dir = os.path.dirname(htsohm.__file__)
    cif_template_path = os.path.join(htsohm_dir, 'simulation', 'forcefield', 'template.cif')
    values = {
            'LatticeConstantA' : round(material.structure.lattice_constant_a, 4),
            'LatticeConstantB' : round(material.structure.lattice_constant_b, 4),
            'LatticeConstantC' : round(material.structure.lattice_constant_c, 4)}
    template_file = open(cif_template_path)
    template = Template( template_file.read() )
    cif_template = template.substitute(values)

    file_name = os.path.join(simulation_path, '%s.cif' % material.uuid)
    with open(file_name, "w") as cif_file:
        cif_file.write(cif_template)
        for atom_site in material.structure.atom_sites:
            cif_file.write(
            "{0:5} C {1:4f} {2:4f} {3:4f}\n".format(
                atom_site.chemical_id,
                round(atom_site.x_frac, 4),
                round(atom_site.y_frac, 4),
                round(atom_site.z_frac, 4)))

def write_mixing_rules(material, simulation_path):
    """Writes .def file for forcefield information.

    """
    htsohm_dir = os.path.dirname(htsohm.__file__)
    mixing_template_path = os.path.join(htsohm_dir, 'simulation', 'forcefield', 'force_field_mixing_rules.def')
    values = {'NumberOfDefinedInteractions' : len(material.structure.lennard_jones) + 10}
    template_file = open(mixing_template_path)
    template = Template( template_file.read() )
    mixing_template = template.substitute(values)

    file_name = os.path.join(simulation_path, 'force_field_mixing_rules.def')
    with open(file_name, "w") as mixing_file:
        mixing_file.write(mixing_template)
        for lennard_jones in material.structure.lennard_jones:
            mixing_file.write(
                "{0:12} lennard-jones {1:8f} {2:8f}\n".format(
                    lennard_jones.chemical_id,
                    round(lennard_jones.epsilon, 4),
                    round(lennard_jones.sigma, 4)
                )
            )
        mixing_file.write(
            "# general mixing rule for Lennard-Jones\n" +
            "Lorentz-Berthelot")

def write_pseudo_atoms(material, simulation_path):
    htsohm_dir = os.path.dirname(htsohm.__file__)
    pseudo_template_path = os.path.join(htsohm_dir, 'simulation', 'forcefield', 'pseudo_atoms.def')
    values = {'NumberOfPseudoAtoms' : len(material.structure.lennard_jones) + 10}
    template_file = open(pseudo_template_path)
    template = Template( template_file.read() )
    pseudo_template = template.substitute(values)

    file_name = os.path.join(simulation_path, 'pseudo_atoms.def')
    with open(file_name, "w") as pseudo_file:
        pseudo_file.write(pseudo_template)
        for atom_type in material.structure.lennard_jones:
            pseudo_file.write(
                "{} yes C C 0 12.0 {} 0.0 1.0 1.0 0 0 absolute 0\n".format(atom_type.chemical_id, 0.0))

def write_force_field(simulation_path):
    """Writes .def file to overwrite LJ-type interactions.

    Args:
        file_name (str): path to .def-file, for example:
            `$(raspa-dir)/forcefield/(run_id)-(uuid)/force_field.def`

    Writes file within RASPA's library:
        `$(raspa-dir)/forcefield/(run_id)-(uuid)/force_field.def`

    NOTE: NO INTERACTIONS ARE OVERWRITTEN BY DEFAULT.

    """
    htsohm_dir = os.path.dirname(htsohm.__file__)
    force_template_path = os.path.join(htsohm_dir, 'simulation', 'forcefield', 'force_field.def')
    template_file = open(force_template_path)
    template = Template( template_file.read() )
    force_template = template.substitute()

    file_name = os.path.join(simulation_path, 'force_field.def')
    with open(file_name, "w") as force_file:
        force_file.write(force_template)
