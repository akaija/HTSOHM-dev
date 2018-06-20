import sys
import os
import subprocess
import shutil
from datetime import datetime
from uuid import uuid4
from string import Template
from pathlib import Path

from htsohm import config
from htsohm.material_files import write_cif_file, write_mixing_rules
from htsohm.material_files import write_pseudo_atoms, write_force_field
from htsohm.simulation.files import load_and_subs_template
from htsohm.db import VoidFraction

def write_raspa_file(filename, seed, simulation_config):
    """Writes RASPA input file for calculating helium void fraction.

    Args:
        filename (str): path to input file.
        run_id (str): identification string for run.
        material_id (str): seed for material.

    Writes RASPA input-file.

    """
    # Load simulation parameters from config
    values = {
            "NumberOfCycles"         : simulation_config["simulation_cycles"],
            "FrameworkName"          : seed,
            "ExternalTemperature"    : simulation_config["temperature"],
            "MoleculeName"           : simulation_config["adsorbate"]}

    # Load template and replace values
    input_data = load_and_subs_template("input_file_templates/void_fraction.input", values)

    # Write simulation input-file
    with open(filename, "w") as raspa_input_file:
        raspa_input_file.write(input_data)

def parse_output(output_file, material, simulation_config):
    """Parse output file for void fraction data.

    Args:
        output_file (str): path to simulation output file.

    Returns:
        results (dict): average Widom Rosenbluth-weight.

    """
    void_fraction = VoidFraction()
    void_fraction.adsorbate = simulation_config["adsorbate"]
    void_fraction.temperature = simulation_config["temperature"]

    with open(output_file) as origin:
        for line in origin:
            if not "Average Widom Rosenbluth-weight:" in line:
                continue
            void_fraction.void_fraction = float(line.split()[4])
        print("\nVOID FRACTION : {}\n".format(void_fraction.void_fraction))

    material.void_fraction.append(void_fraction)

def run(material, structure, simulation_config):
    """Runs void fraction simulation.

    Args:
        material (Material): material record.

    Returns:
        results (dict): void fraction simulation results.

    """
    # Determine where to write simulation input/output files, create directory
    simulation_directory  = config["simulation_directory"]
    if simulation_directory == "HTSOHM":
        htsohm_dir = config["htsohm_dir"]
        path = os.path.join(htsohm_dir, material.run_id)
    elif simulation_directory == "SCRATCH":
        path = os.environ["SCRATCH"]
    else:
        print("OUTPUT DIRECTORY NOT FOUND.")
    output_dir = os.path.join(path, "output_{}_{}".format(material.seed, uuid4()))
    print("Output directory : {}".format(output_dir))
    os.makedirs(output_dir, exist_ok=True)

    # Write simulation input-files
    # RASPA input-file
    filename = os.path.join(output_dir, "VoidFraction.input")
    write_raspa_file(filename, material.seed, simulation_config)
    # Pseudomaterial cif-file
    write_cif_file(material, structure, output_dir)
    # Lennard-Jones parameters, force_field_mixing_rules.def
    write_mixing_rules(structure, output_dir)
    # Pseudoatom definitions, pseudo_atoms.def (placeholder values)
    write_pseudo_atoms(structure, output_dir)
    # Overwritten interactions, force_field.def (none overwritten by default)
    write_force_field(output_dir)

    # Run simulations
    while True:
        try:
            print("Date             : {}".format(datetime.now().date().isoformat()))
            print("Time             : {}".format(datetime.now().time().isoformat()))
            print("Simulation type  : {}".format(simulation_config["type"]))
            print("Probe            : {}".format(simulation_config["adsorbate"]))
            print("Temperature      : {}".format(simulation_config["temperature"]))
            filename = "output_{}_2.2.2_298.000000_0.data".format(material.seed)
            output_file = os.path.join(output_dir, "Output", "System_0", filename)
            while not Path(output_file).exists():
                process = subprocess.run(["simulate", "./VoidFraction.input"], check=True, cwd=output_dir)

            # Parse output
            parse_output(output_file, material, simulation_config)
            shutil.rmtree(output_dir, ignore_errors=True)
            sys.stdout.flush()
        except (FileNotFoundError, IndexError, KeyError) as err:
            print(err)
            print(err.args)
            continue
        break
