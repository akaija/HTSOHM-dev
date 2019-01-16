import os

def write_mol_file(material, structure, simulation_path):
    """Writes .mol file for structural information."""
    file_name = os.path.join(simulation_path, "{}.mol".format(material.seed))
    with open(file_name, "w") as mol_file:
        mol_file.write(
                " Molecule_name: {}\n".format(material.seed) +
                "\n" +
                "  Coord_Info: Listed Cartesian None\n" +
                "        {}\n".format(structure.n()))
        for i in range(structure.n()):
            a = structure.atom_sites[i]
            mol_file.write(
                    "{:6} {:10.4f} {:10.4f} {:10.4f}  {:5} {:10.8f}  0  0\n".format(
                        i + 1, round(structure.atom_sites[i].x * structure.lattice_constants.a, 4),
                        round(structure.atom_sites[i].y * structure.lattice_constants.b, 4),
                        round(structure.atom_sites[i].z * structure.lattice_constants.c, 4),
                        structure.atom_sites[i].chemical_id, round(structure.atom_sites[i].q, 8)))
        mol_file.write(
                "\n" +
                "\n" +
                "\n" +
                "  Fundcell_Info: Listed\n" +
                "        {:10.4f}       {:10.4f}       {:10.4f}\n".format(
                    round(structure.lattice_constants.a, 4),
                    round(structure.lattice_constants.b, 4),
                    round(structure.lattice_constants.c, 4)) +
                "           90.0000          90.0000          90.0000\n" +
                "           0.00000          0.00000          0.00000\n" +
                "        {:10.4f}       {:10.4f}       {:10.4f}\n".format(
                    round(structure.lattice_constants.a, 4),
                    round(structure.lattice_constants.b, 4),
                    round(structure.lattice_constants.c, 4)) +
                "\n")

def write_mixing_rules(structure, simulation_path):
    """Writes .def file for forcefield information.

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
            "{}\n".format(len(structure.atom_types) + 10) +
            "# type interaction, parameters.    " +
            "IMPORTANT: define shortest matches first, so" +
            " that more specific ones overwrites these\n"
        )
        for lj in structure.atom_types:
            mixing_rules_file.write(
                "{0:12} lennard-jones {1:8f} {2:8f}\n".format(lj.chemical_id,
                    round(lj.epsilon, 4), round(lj.sigma, 4)))
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

def write_pseudo_atoms(structure, simulation_path):
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
            "%s\n" % (len(structure.atom_types) + 10) +
            "#type  print   as  chem    oxidation   mass    charge  polarization    B-factor    radii   " +
                 "connectivity  anisotropic anisotrop-type  tinker-type\n")
        for a in structure.atom_types:
            pseudo_atoms_file.write(
                "{0:7}  yes  C   C   0   12.0       0.0  0.0  1.0  1.0    0  0  absolute  0\n".format(
                    a.chemical_id))
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
