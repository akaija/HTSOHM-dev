import sys
import os
import subprocess
import shutil
from datetime import datetime
from uuid import uuid4
from string import Template

import htsohm
from htsohm.simulation.forcefield.files import write_cif_file, write_mixing_rules
from htsohm.simulation.forcefield.files import write_force_field, write_pseudo_atoms
from htsohm.simulation.calculate_bin import calc_bin

def write_raspa_file(filename, uuid, params):
    """Writes RASPA input file for calculating helium void fraction.

    Args:
        filename (str): path to input file.
        run_id (str): identification string for run.
        material_id (str): uuid for material.

    Writes RASPA input-file.

    """
    # Load simulation parameters from config
    values = {
            'NumberOfCycles'                : params['simulation_cycles'],
            'FrameworkName'                 : uuid}

    # Load template and replace values
    htsohm_dir = os.path.dirname(htsohm.__file__)
    input_template_path = os.path.join(htsohm_dir, 'simulation', 'helium_void_fraction.input')
    input_file = open(input_template_path)
    template = Template( input_file.read() )
    input_data = template.substitute(values)

    print(input_data)

    # Write simulation input-file
    with open(filename, "w") as raspa_input_file:
        raspa_input_file.write(input_data)

def parse_output(output_file, params):
    """Parse output file for void fraction data.

    Args:
        output_file (str): path to simulation output file.

    Returns:
        results (dict): average Widom Rosenbluth-weight.

    """
    results = {}
    with open(output_file) as origin:
        for line in origin:
            if not "Average Widom Rosenbluth-weight:" in line:
                continue
            results['vf_helium_void_fraction'] = float(line.split()[4])
        print("\nVOID FRACTION :   %s\n" % (results['vf_helium_void_fraction']))

    # calculate bin
    results['void_fraction_bin'] = calc_bin(
                results['vf_helium_void_fraction'],
                *params['limits'],
                params['bins'])

    return results

def run(material, simulation_config):
    """Runs void fraction simulation.

    Args:
        material (Material): material record.

    Returns:
        results (dict): void fraction simulation results.

    """
    # Determine where to write simulation input/output files, create directory
    simulation_directory  = simulation_config['directory']
    if simulation_directory == 'HTSOHM':
        htsohm_dir = os.path.dirname(os.path.dirname(htsohm.__file__))
        path = os.path.join(htsohm_dir, material.run_id)
    elif simulation_directory == 'SCRATCH':
        path = os.environ['SCRATCH']
    else:
        print('OUTPUT DIRECTORY NOT FOUND.')
    output_dir = os.path.join(path, 'output_%s_%s' % (material.uuid, uuid4()))
    print("Output directory :\t%s" % output_dir)
    os.makedirs(output_dir, exist_ok=True)

    # Write simulation input-files
    params = simulation_config['helium_void_fraction']
    # RASPA input-file
    filename = os.path.join(output_dir, "VoidFraction.input")
    write_raspa_file(filename, material.uuid, params)
    # Pseudomaterial cif-file
    write_cif_file(material, output_dir)
    # Lennard-Jones parameters, force_field_mixing_rules.def
    write_mixing_rules(material, output_dir)
    # Pseudoatom definitions, pseudo_atoms.def
    write_pseudo_atoms(material, output_dir)
    # Overwritten interactions, force_field.def (none overwritten by default)
    write_force_field(output_dir)

    # Run simulations
    while True:
        try:
            print("Date :\t%s" % datetime.now().date().isoformat())
            print("Time :\t%s" % datetime.now().time().isoformat())
            print("Calculating void fraction of %s..." % (material.uuid))
            subprocess.run(['simulate', './VoidFraction.input'], check=True, cwd=output_dir)
            filename = "output_%s_1.1.1_298.000000_0.data" % (material.uuid)
            output_file = os.path.join(output_dir, 'Output', 'System_0', filename)

            # Parse output
            results = parse_output(output_file, params)
            shutil.rmtree(output_dir, ignore_errors=True)
            sys.stdout.flush()
        except (FileNotFoundError, IndexError, KeyError) as err:
            print(err)
            print(err.args)
            continue
        break

    return results
